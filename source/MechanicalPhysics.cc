/* SPDX-FileCopyrightText: Copyright (c) 2022 - 2026, the adamantine authors.
 * SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
 */

#include <MechanicalPhysics.hh>
#include <instantiation.hh>

#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/numerics/vector_tools.h>

#ifdef ADAMANTINE_WITH_CALIPER
#include <caliper/cali.h>
#endif

// FOR DEBUGGING
#include <iostream>
#include <algorithm>
#include <cmath>
#include <map>
#include <exception>
#include <fstream>
#include <limits>
#include <deal.II/lac/precondition.h>
// --

namespace
{
template <typename VectorType>
double
project_out_nullspace(VectorType &vector,
                      std::vector<VectorType> const &nullspace_modes)
{
  double removed_norm_squared = 0.0;

  for (auto const &mode : nullspace_modes)
  {
    double const coefficient = mode * vector;
    removed_norm_squared += coefficient * coefficient;
    vector.add(-coefficient, mode);
  }

  return std::sqrt(removed_norm_squared);
}

template <typename MatrixType, typename VectorType>
class NullspaceProjectedOperator
{
public:
  NullspaceProjectedOperator(
      MatrixType const &matrix,
      std::vector<VectorType> const &nullspace_modes,
      dealii::IndexSet const &locally_owned_dofs,
      MPI_Comm const communicator)
      : _matrix(matrix), _nullspace_modes(nullspace_modes),
        _scratch(locally_owned_dofs, communicator)
  {
  }

  void vmult(VectorType &dst, VectorType const &src) const
  {
    _scratch = src;
    project_out_nullspace(_scratch, _nullspace_modes);

    _matrix.vmult(dst, _scratch);
    project_out_nullspace(dst, _nullspace_modes);
  }

private:
  MatrixType const &_matrix;
  std::vector<VectorType> const &_nullspace_modes;
  mutable VectorType _scratch;
};

template <typename PreconditionerType, typename VectorType>
class NullspaceProjectedPreconditioner
{
public:
  NullspaceProjectedPreconditioner(
      PreconditionerType const &preconditioner,
      std::vector<VectorType> const &nullspace_modes,
      dealii::IndexSet const &locally_owned_dofs,
      MPI_Comm const communicator)
      : _preconditioner(preconditioner), _nullspace_modes(nullspace_modes),
        _scratch(locally_owned_dofs, communicator)
  {
  }

  void vmult(VectorType &dst, VectorType const &src) const
  {
    _scratch = src;
    project_out_nullspace(_scratch, _nullspace_modes);

    _preconditioner.vmult(dst, _scratch);
    project_out_nullspace(dst, _nullspace_modes);
  }

private:
  PreconditionerType const &_preconditioner;
  std::vector<VectorType> const &_nullspace_modes;
  mutable VectorType _scratch;
};
} // namespace

namespace adamantine
{
template <int dim, int n_materials, int p_order, typename MaterialStates,
          typename MemorySpaceType>
MechanicalPhysics<dim, n_materials, p_order, MaterialStates, MemorySpaceType>::
    MechanicalPhysics(
        MPI_Comm const &communicator, unsigned int const fe_degree,
        Geometry<dim> &geometry, Boundary const &boundary,
        MaterialProperty<dim, n_materials, p_order, MaterialStates,
                         MemorySpaceType> &material_properties,
        std::vector<double> const &reference_temperatures)
    : _geometry(geometry), _boundary(boundary),
      _material_properties(material_properties),
      _dof_handler(_geometry.get_triangulation()),
      _reference_temperatures(reference_temperatures),
      _solution_transfer(_dof_handler),
      _closest_quad_point_adaptation(dealii::QGauss<dim>(fe_degree + 1)),
      _cell_data_transfer(
          dynamic_cast<const dealii::parallel::distributed::Triangulation<dim>
                           &>(_dof_handler.get_triangulation()),
          /* transfer_variable_size_data */ false,
          [&](const typename dealii::Triangulation<dim>::cell_iterator &parent,
              const std::vector<std::vector<double>> &parent_values)
          {
            return _closest_quad_point_adaptation.coarse_to_fine(parent,
                                                                 parent_values);
          },
          [&](const typename dealii::Triangulation<dim>::cell_iterator &parent,
              const std::vector<std::vector<std::vector<double>>> &child_values)
          {
            return _closest_quad_point_adaptation.fine_to_coarse(parent,
                                                                 child_values);
          })
{
  // Create the FECollection
  _fe_collection.push_back(
      dealii::FESystem<dim>(dealii::FE_Q<dim>(fe_degree) ^ dim));
  _fe_collection.push_back(
      dealii::FESystem<dim>(dealii::FE_Nothing<dim>() ^ dim));

  // Create the QCollection
  _q_collection.push_back(dealii::QGauss<dim>(fe_degree + 1));
  _q_collection.push_back(dealii::QGauss<dim>(1));

  // Solve the mechanical problem only on the part of the domain that has solid
  // material.
  unsigned int n_active_cells =
      _dof_handler.get_triangulation().n_active_cells();
  for (auto const &cell :
       dealii::filter_iterators(_dof_handler.active_cell_iterators(),
                                dealii::IteratorFilters::LocallyOwnedCell()))
  {
    if (_material_properties.get_state_ratio(
            cell, MaterialStates::State::solid) > 0.99)
    {
      cell->set_active_fe_index(0);
    }
    else
    {
      cell->set_active_fe_index(1);
    }
  }

  // Create the mechanical operator
  _mechanical_operator =
      std::make_unique<MechanicalOperator<dim, n_materials, p_order,
                                          MaterialStates, MemorySpaceType>>(
          communicator, _material_properties, reference_temperatures);

  // Create the data used to compute the stress tensor
  unsigned int const n_quad_pts = _q_collection.max_n_quadrature_points();
  _plastic_internal_variable.reserve(n_active_cells);

  for (auto const &cell : _dof_handler.active_cell_iterators())
  {
    if (cell->is_locally_owned())
    {
      auto elastic_limit = _material_properties.get_mechanical_property(
          cell, StateProperty::elastic_limit);
      _plastic_internal_variable.emplace_back(
          std::vector<double>(n_quad_pts, elastic_limit));
    }
    else
    {
      _plastic_internal_variable.emplace_back(std::vector<double>(
          n_quad_pts, std::numeric_limits<double>::signaling_NaN()));
    }
  }
  _stress.resize(n_active_cells,
                 std::vector<dealii::SymmetricTensor<2, dim>>(n_quad_pts));
  _back_stress.resize(n_active_cells,
                      std::vector<dealii::SymmetricTensor<2, dim>>(n_quad_pts));
  _thermal_stress.resize(n_active_cells, std::vector<double>(n_quad_pts));
}

template <int dim, int n_materials, int p_order, typename MaterialStates,
          typename MemorySpaceType>
void MechanicalPhysics<dim, n_materials, p_order, MaterialStates,
                       MemorySpaceType>::
    setup_dofs(std::vector<std::shared_ptr<BodyForce<dim>>> const &body_forces)
{
  _dof_handler.distribute_dofs(_fe_collection);
  dealii::IndexSet locally_relevant_dofs =
      dealii::DoFTools::extract_locally_relevant_dofs(_dof_handler);
  dealii::IndexSet locally_owned_dofs = _dof_handler.locally_owned_dofs();
  _affine_constraints.reinit(locally_owned_dofs, locally_relevant_dofs);
  dealii::DoFTools::make_hanging_node_constraints(_dof_handler,
                                                  _affine_constraints);

  std::map<dealii::types::boundary_id, const dealii::Function<dim> *>
      boundary_function_map;
  dealii::Functions::ZeroFunction<dim> zero_function(dim);
  auto boundary_ids = _boundary.get_boundary_ids(BoundaryType::clamped);
  
  for (auto id : boundary_ids)
  {
    boundary_function_map[id] = &zero_function;
  }
  dealii::VectorTools::interpolate_boundary_values(
      _dof_handler, boundary_function_map, _affine_constraints);
  _affine_constraints.close();

  _mechanical_operator->reinit(_dof_handler, _affine_constraints, _q_collection,
                               body_forces);
}

template <int dim, int n_materials, int p_order, typename MaterialStates,
          typename MemorySpaceType>
void MechanicalPhysics<dim, n_materials, p_order, MaterialStates,
                       MemorySpaceType>::
    update_rhs(std::vector<std::shared_ptr<BodyForce<dim>>> const &body_forces)
{
  _mechanical_operator->assemble_rhs(body_forces);
}

template <int dim, int n_materials, int p_order, typename MaterialStates,
          typename MemorySpaceType>
void MechanicalPhysics<dim, n_materials, p_order, MaterialStates,
                       MemorySpaceType>::
    update_rhs(
        dealii::DoFHandler<dim> const &thermal_dof_handler,
        dealii::LA::distributed::Vector<double, dealii::MemorySpace::Host> const
            &temperature,
        std::vector<bool> const &has_melted,
        std::vector<std::shared_ptr<BodyForce<dim>>> const &body_forces)
{
  _thermal_dof_handler = &thermal_dof_handler;
  _temperature = temperature;
  _has_melted = has_melted;
  _mechanical_operator->update_temperature(thermal_dof_handler, temperature,
                                           has_melted);
  _mechanical_operator->assemble_rhs(body_forces);
}

template <int dim, int n_materials, int p_order, typename MaterialStates,
          typename MemorySpaceType>
void MechanicalPhysics<dim, n_materials, p_order, MaterialStates,
                       MemorySpaceType>::prepare_transfer_mpi()
{
  _old_displacement.update_ghost_values();
  _solution_transfer.prepare_for_coarsening_and_refinement(_old_displacement);

  _data_to_transfer.clear();
  unsigned int const n_quad_pts = _q_collection.max_n_quadrature_points();
  unsigned int const n_doubles_per_quad_scalar = 2;
  unsigned int const n_doubles_per_quad_stress =
      dealii::SymmetricTensor<2, dim>::n_independent_components;

  unsigned int const n_doubles_per_quad =
      n_doubles_per_quad_scalar + n_doubles_per_quad_stress * 2;
  std::vector<std::vector<double>> dummy_cell_data(
      n_quad_pts, std::vector<double>(n_doubles_per_quad,
                                      std::numeric_limits<double>::infinity()));
  std::vector<std::vector<double>> cell_data = dummy_cell_data;
  unsigned int cell_id = 0;
  for (auto const &cell : _dof_handler.active_cell_iterators())
  {
    if (cell->is_locally_owned())
    {
      unsigned int const stress_offset = n_doubles_per_quad_scalar;
      unsigned int const back_stress_offset =
          n_doubles_per_quad_scalar + n_doubles_per_quad_stress;

      for (unsigned int quad = 0; quad < n_quad_pts; ++quad)
      {
        std::vector<double> &cell_data_quad = cell_data[quad];
        cell_data_quad[0] = _plastic_internal_variable[cell_id][quad];
        cell_data_quad[1] = _thermal_stress[cell_id][quad];
        for (unsigned int i = 0; i < n_doubles_per_quad_stress; ++i)
        {
          cell_data_quad[stress_offset + i] =
              _stress[cell_id][quad].access_raw_entry(i);
          cell_data_quad[back_stress_offset + i] =
              _back_stress[cell_id][quad].access_raw_entry(i);
        }
      }
      _data_to_transfer.push_back(cell_data);
    }
    else
    {
      _data_to_transfer.push_back(dummy_cell_data);
    }
    ++cell_id;
  }
  _cell_data_transfer.prepare_for_coarsening_and_refinement(_data_to_transfer);
}

template <int dim, int n_materials, int p_order, typename MaterialStates,
          typename MemorySpaceType>
void MechanicalPhysics<dim, n_materials, p_order, MaterialStates,
                       MemorySpaceType>::complete_transfer_mpi()
{
  _dof_handler.distribute_dofs(_fe_collection);

  const dealii::IndexSet locally_relevant_dofs =
      dealii::DoFTools::extract_locally_relevant_dofs(_dof_handler);
  _old_displacement.reinit(_dof_handler.locally_owned_dofs(),
                           locally_relevant_dofs,
#if DEAL_II_VERSION_GTE(9, 7, 0)
                           _dof_handler.get_mpi_communicator()
#else
                           _dof_handler.get_communicator()
#endif
  );
  _solution_transfer.interpolate(_old_displacement);

  auto n_active_cells = _dof_handler.get_triangulation().n_active_cells();
  unsigned int const n_quad_pts = _q_collection.max_n_quadrature_points();

  _plastic_internal_variable.resize(n_active_cells,
                                    std::vector<double>(n_quad_pts));
  _thermal_stress.resize(n_active_cells, std::vector<double>(n_quad_pts));
  _stress.resize(n_active_cells,
                 std::vector<dealii::SymmetricTensor<2, dim>>(n_quad_pts));
  _back_stress.resize(n_active_cells,
                      std::vector<dealii::SymmetricTensor<2, dim>>(n_quad_pts));

  unsigned int const n_doubles_per_quad_scalar = 2;
  unsigned int const n_doubles_per_quad_stress =
      dealii::SymmetricTensor<2, dim>::n_independent_components;

  unsigned int const n_doubles_per_quad =
      n_doubles_per_quad_scalar + n_doubles_per_quad_stress * 2;
  std::vector<std::vector<std::vector<double>>> data_to_unpack(
      n_active_cells, std::vector<std::vector<double>>(
                          n_quad_pts, std::vector<double>(n_doubles_per_quad)));
  _cell_data_transfer.unpack(data_to_unpack);

  unsigned int cell_id = 0;
  for (auto const &cell : _dof_handler.active_cell_iterators())
  {
    if (cell->is_locally_owned())
    {
      unsigned int const stress_offset = n_doubles_per_quad_scalar;
      unsigned int const back_stress_offset =
          n_doubles_per_quad_scalar + n_doubles_per_quad_stress;

      for (unsigned int quad = 0; quad < n_quad_pts; ++quad)
      {
        _plastic_internal_variable[cell_id][quad] =
            data_to_unpack[cell_id][quad][0];
        _thermal_stress[cell_id][quad] = data_to_unpack[cell_id][quad][1];
        for (unsigned int i = 0; i < n_doubles_per_quad_stress; ++i)
        {
          _stress[cell_id][quad].access_raw_entry(i) =
              data_to_unpack[cell_id][quad][stress_offset + i];
          _back_stress[cell_id][quad].access_raw_entry(i) =
              data_to_unpack[cell_id][quad][back_stress_offset + i];
        }
      }
    }
    ++cell_id;
  }
}

template <int dim, int n_materials, int p_order, typename MaterialStates,
          typename MemorySpaceType>
void MechanicalPhysics<dim, n_materials, p_order, MaterialStates,
                       MemorySpaceType>::
    setup_dofs(
        dealii::DoFHandler<dim> const &thermal_dof_handler,
        dealii::LA::distributed::Vector<double, dealii::MemorySpace::Host> const
            &temperature,
        std::vector<bool> const &has_melted, bool rebuild_matrix,
        std::vector<std::shared_ptr<BodyForce<dim>>> const &body_forces)
{
  _thermal_dof_handler = &thermal_dof_handler;
  _temperature = temperature;
  _has_melted = has_melted;
  _mechanical_operator->update_temperature(thermal_dof_handler, temperature,
                                           has_melted);
  // Update the active fe indices, the plastic variables, and the displacement.
  unsigned int const n_quad_pts = _q_collection.max_n_quadrature_points();
  unsigned int cell_id = 0;
  std::vector<std::vector<double>> saved_old_displacement;
  std::vector<std::vector<double>> tmp_plastic_internal_variable;
  std::vector<std::vector<double>> tmp_thermal_stress;
  std::vector<std::vector<dealii::SymmetricTensor<2, dim>>> tmp_stress;
  std::vector<std::vector<dealii::SymmetricTensor<2, dim>>> tmp_back_stress;
  // The number of cells to activate/deactive should be small, so we can
  // already reserve the memory.
  unsigned int const n_dofs_per_cell = _fe_collection.max_dofs_per_cell();
  unsigned int const n_old_active_cells = _plastic_internal_variable.size();
  std::vector<dealii::types::global_dof_index> global_dof_indices(
      n_dofs_per_cell);
  tmp_plastic_internal_variable.reserve(n_old_active_cells);
  tmp_thermal_stress.reserve(n_old_active_cells);
  tmp_stress.reserve(n_old_active_cells);
  tmp_back_stress.reserve(_back_stress.size());
  // First we save _old_displacement if it exists
  if (_old_displacement.size())
  {
    _old_displacement.update_ghost_values();

    std::vector<double> cell_values(n_dofs_per_cell);
    saved_old_displacement.reserve(n_old_active_cells);
    for (auto const &cell : _dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
      {
        auto fe_index = cell->active_fe_index();
        if (fe_index == 0)
        {
          // The cell contains solid material, we need to save the displacement
          cell->get_dof_indices(global_dof_indices);
          for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
          {
            cell_values[i] = _old_displacement[global_dof_indices[i]];
          }
        }
        else
        {
          // The cell does not contain material or it is liquid. The
          // displacement is ignored.
          cell_values.assign(n_dofs_per_cell, 0.);
        }
        saved_old_displacement.push_back(cell_values);
      }
      else
      {
        saved_old_displacement.push_back(std::vector<double>(n_dofs_per_cell));
      }
    }
  }

  // Now we can update the fe indices and the plastic variables.
  for (auto const &cell : _dof_handler.active_cell_iterators())
  {
    if (cell->is_locally_owned())
    {
      auto current_fe_index = cell->active_fe_index();
      if (_material_properties.get_state_ratio(
              cell, MaterialStates::State::solid) > 0.99)
      {
        // Only enable the cell if it is also enabled for the thermal simulation
        // Get the thermal DoFHandler cell iterator
        dealii::DoFCellAccessor<dim, dim, false> thermal_cell(
            &(_dof_handler.get_triangulation()), cell->level(), cell->index(),
            &thermal_dof_handler);
        auto updated_fe_index = thermal_cell.active_fe_index();
        if (current_fe_index == updated_fe_index)
        {
          // The cells is unchanged, we just copy the plastic variables as-is.
          tmp_plastic_internal_variable.push_back(
              _plastic_internal_variable[cell_id]);
          tmp_thermal_stress.push_back(_thermal_stress[cell_id]);
          tmp_stress.push_back(_stress[cell_id]);
          tmp_back_stress.push_back(_back_stress[cell_id]);
        }
        else
        {
          // The cell has solidified or material has been added. The new cells
          // are initialized with default values.
          auto elastic_limit = _material_properties.get_mechanical_property(
              cell, StateProperty::elastic_limit);
          tmp_plastic_internal_variable.push_back(
              std::vector<double>(n_quad_pts, elastic_limit));
          tmp_thermal_stress.push_back(std::vector<double>(n_quad_pts));
          tmp_stress.push_back(
              std::vector<dealii::SymmetricTensor<2, dim>>(n_quad_pts));
          tmp_back_stress.push_back(
              std::vector<dealii::SymmetricTensor<2, dim>>(n_quad_pts));

          cell->set_active_fe_index(updated_fe_index);
          rebuild_matrix = true;
        }
      }
      else
      {
        if (current_fe_index == 0)
        {
          rebuild_matrix = true;
        }

        // The cell is liquid. We don't need to save the plastic variables.
        cell->set_active_fe_index(1);
        tmp_plastic_internal_variable.push_back(std::vector<double>(
            n_quad_pts, std::numeric_limits<double>::signaling_NaN()));
        tmp_thermal_stress.push_back(std::vector<double>(n_quad_pts));
        tmp_stress.push_back(
            std::vector<dealii::SymmetricTensor<2, dim>>(n_quad_pts));
        tmp_back_stress.push_back(
            std::vector<dealii::SymmetricTensor<2, dim>>(n_quad_pts));
      }
    }
    else
    {
      tmp_plastic_internal_variable.push_back(std::vector<double>(
          n_quad_pts, std::numeric_limits<double>::signaling_NaN()));
      tmp_thermal_stress.push_back(std::vector<double>(n_quad_pts));
      tmp_stress.push_back(
          std::vector<dealii::SymmetricTensor<2, dim>>(n_quad_pts));
      tmp_back_stress.push_back(
          std::vector<dealii::SymmetricTensor<2, dim>>(n_quad_pts));
    }
    ++cell_id;
  }

  // Check if we need to rebuild the matrix
  rebuild_matrix =
      dealii::Utilities::MPI::logical_or(rebuild_matrix,
#if DEAL_II_VERSION_GTE(9, 7, 0)
                                         _dof_handler.get_mpi_communicator()
#else
                                         _dof_handler.get_communicator()
#endif
      );
  // If we do not need to rebuild the matrix. Update the rhs and exit.
  if (!rebuild_matrix)
  {
    update_rhs(body_forces);
    return;
  }

  _plastic_internal_variable.swap(tmp_plastic_internal_variable);
  _thermal_stress.swap(tmp_thermal_stress);
  _stress.swap(tmp_stress);
  _back_stress.swap(tmp_back_stress);

  setup_dofs(body_forces);

  // Update _old_displacement if necessary
  const dealii::IndexSet locally_relevant_dofs =
      dealii::DoFTools::extract_locally_relevant_dofs(_dof_handler);
  const dealii::IndexSet locally_owned_dofs = _dof_handler.locally_owned_dofs();
  _old_displacement.reinit(locally_owned_dofs, locally_relevant_dofs,
#if DEAL_II_VERSION_GTE(9, 7, 0)
                           _dof_handler.get_mpi_communicator()
#else
                           _dof_handler.get_communicator()
#endif
  );

  if (saved_old_displacement.size())
  {
    cell_id = 0;
    for (auto const &cell : _dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned())
      {
        auto fe_index = cell->active_fe_index();
        if (fe_index == 0)
        {
          cell->get_dof_indices(global_dof_indices);
          for (unsigned int i = 0; i < n_dofs_per_cell; ++i)
          {
            if (locally_owned_dofs.is_element(global_dof_indices[i]))
              _old_displacement[global_dof_indices[i]] =
                  saved_old_displacement[cell_id][i];
          }
        }
      }
      ++cell_id;
    }
    _old_displacement.compress(dealii::VectorOperation::insert);
  }
}

template <int dim, int n_materials, int p_order, typename MaterialStates,
          typename MemorySpaceType>
void
MechanicalPhysics<dim, n_materials, p_order, MaterialStates,
                  MemorySpaceType>::detect_floating_rigid_body_modes()
{
  _floating_rigid_body_modes.clear();

#if DEAL_II_VERSION_GTE(9, 7, 0)
  MPI_Comm const mechanics_communicator =
      _dof_handler.get_mpi_communicator();
#else
  MPI_Comm const mechanics_communicator =
      _dof_handler.get_communicator();
#endif

  unsigned int const n_mpi_processes =
      dealii::Utilities::MPI::n_mpi_processes(mechanics_communicator);

  // Recompute hanging-node constraints separately. At this point,
  // _affine_constraints also contains the clamped boundary constraints,
  // so _affine_constraints.n_constraints() cannot be used to determine
  // whether hanging-node constraints are present.
  dealii::IndexSet const locally_owned_dofs =
      _dof_handler.locally_owned_dofs();

  dealii::IndexSet const locally_relevant_dofs =
      dealii::DoFTools::extract_locally_relevant_dofs(_dof_handler);

  dealii::AffineConstraints<double> hanging_constraints;
  hanging_constraints.reinit(locally_owned_dofs, locally_relevant_dofs);

  dealii::DoFTools::make_hanging_node_constraints(
      _dof_handler, hanging_constraints);

  auto const n_hanging_constraints =
      hanging_constraints.n_constraints();

  auto const boundary_ids =
      _boundary.get_boundary_ids(BoundaryType::clamped);

  if (n_mpi_processes == 1 && n_hanging_constraints == 0)
  {
    using CellIterator =
        typename dealii::DoFHandler<dim>::active_cell_iterator;

    std::vector<CellIterator> mechanical_cells;

    for (auto const &cell : _dof_handler.active_cell_iterators())
    {
      if (cell->is_locally_owned() && cell->active_fe_index() == 0)
        mechanical_cells.push_back(cell);
    }

    unsigned int const n_mechanical_cells =
        mechanical_cells.size();

    // Union-find data structure used to identify face-connected
    // mechanically active components.
    std::vector<unsigned int> parent(n_mechanical_cells);

    for (unsigned int i = 0; i < n_mechanical_cells; ++i)
      parent[i] = i;

    auto find_root = [&parent](unsigned int i)
    {
      while (parent[i] != i)
      {
        parent[i] = parent[parent[i]];
        i = parent[i];
      }

      return i;
    };

    auto unite = [&parent, &find_root](unsigned int a,
                                       unsigned int b)
    {
      a = find_root(a);
      b = find_root(b);

      if (a != b)
        parent[b] = a;
    };

    // Keep track of which mechanical cell first owns each interior face.
    // When another active mechanical cell sees the same face, the two cells
    // belong to the same connected component.
    std::map<unsigned int, unsigned int> face_owner;

    // Track whether each mechanical cell touches one of the clamped
    // geometric boundaries.
    std::vector<bool> cell_touches_clamp(
        n_mechanical_cells, false);

    for (unsigned int i = 0; i < n_mechanical_cells; ++i)
    {
      auto const &cell = mechanical_cells[i];

      for (unsigned int f = 0;
           f < dealii::GeometryInfo<dim>::faces_per_cell;
           ++f)
      {
        auto const face = cell->face(f);

        if (face->at_boundary())
        {
          auto const boundary_id = face->boundary_id();

          if (std::find(boundary_ids.begin(),
                        boundary_ids.end(),
                        boundary_id) != boundary_ids.end())
          {
            cell_touches_clamp[i] = true;
          }
        }
        else
        {
          unsigned int const face_index = face->index();

          auto const result =
              face_owner.emplace(face_index, i);

          if (!result.second)
            unite(i, result.first->second);
        }
      }
    }

    struct ComponentInfo
    {
      std::vector<unsigned int> cell_indices;
      bool touches_clamped_boundary = false;
    };

    // Collapse the cell-level union-find information into connected
    // mechanical components.
    std::map<unsigned int, ComponentInfo> components;

    for (unsigned int i = 0; i < n_mechanical_cells; ++i)
    {
      unsigned int const root = find_root(i);

      components[root].cell_indices.push_back(i);

      components[root].touches_clamped_boundary =
          components[root].touches_clamped_boundary ||
          cell_touches_clamp[i];
    }

    // deal.II 9.6.x does not provide
    // DoFTools::extract_rigid_body_modes(), so construct the rigid-body
    // displacement fields explicitly.
    //
    // DoFTools::map_dofs_to_support_points() gives the real-space
    // coordinate associated with each displacement DoF. For r measured
    // from the floating-component center, the rigid modes are translations
    // and omega x r rotations.
    dealii::MappingQ1<dim> mapping;

    auto const support_points =
        dealii::DoFTools::map_dofs_to_support_points(
            mapping, _dof_handler);

    std::vector<dealii::types::global_dof_index>
        local_dof_indices(
            _dof_handler.get_fe_collection().max_dofs_per_cell());

    unsigned int floating_component_number = 0;

    for (auto const &[root, component] : components)
    {
      (void)root;

      // This component is mechanically supported by a clamped boundary,
      // so it does not contribute an independent rigid-body nullspace.
      if (component.touches_clamped_boundary)
        continue;

      // Record which global DoFs belong to this floating component and
      // which displacement component (x/y/z) each DoF represents.
      //
      // The active mechanical element is FESystem(FE_Q ^ dim), so every
      // shape function is primitive and system_to_component_index() is
      // well defined.
      std::vector<bool> component_dofs(
          _dof_handler.n_dofs(), false);

      std::vector<unsigned int> dof_components(
          _dof_handler.n_dofs(), dim);

      for (auto const cell_index : component.cell_indices)
      {
        auto const &cell = mechanical_cells[cell_index];
        auto const &fe = cell->get_fe();

        unsigned int const dofs_per_cell =
            fe.n_dofs_per_cell();

        cell->get_dof_indices(local_dof_indices);

        for (unsigned int i = 0;
             i < dofs_per_cell;
             ++i)
        {
          auto const dof = local_dof_indices[i];

          component_dofs[dof] = true;

          dof_components[dof] =
              fe.system_to_component_index(i).first;
        }
      }

      // Express rotations about the component center rather than the
      // global origin. A change of rotation center only adds a translation,
      // so the nullspace span is unchanged. Centering improves numerical
      // scaling before Gram-Schmidt orthonormalization.
      dealii::Point<dim> component_center;

      unsigned int n_component_dofs_with_support_point = 0;

      for (auto const dof : locally_owned_dofs)
      {
        if (!component_dofs[dof])
          continue;

        auto const point = support_points.find(dof);

        if (point == support_points.end())
          continue;

        for (unsigned int d = 0; d < dim; ++d)
          component_center[d] += point->second[d];

        ++n_component_dofs_with_support_point;
      }

      if (n_component_dofs_with_support_point == 0)
      {
        std::cout
            << "WARNING: floating mechanical component "
            << floating_component_number
            << " has no mapped support points; rigid-mode "
               "projection is skipped for this component."
            << std::endl;

        ++floating_component_number;
        continue;
      }

      for (unsigned int d = 0; d < dim; ++d)
      {
        component_center[d] /=
            static_cast<double>(
                n_component_dofs_with_support_point);
      }

      // Number of rigid motions:
      //
      // 1D: 1 translation
      // 2D: 2 translations + 1 rotation = 3
      // 3D: 3 translations + 3 rotations = 6
      unsigned int const n_rigid_body_modes =
          dim + dim * (dim - 1) / 2;

      unsigned int modes_added_for_component = 0;

      for (unsigned int rigid_mode = 0;
           rigid_mode < n_rigid_body_modes;
           ++rigid_mode)
      {
        dealii::LA::distributed::Vector<
            double, dealii::MemorySpace::Host>
            mode(locally_owned_dofs,
                 locally_relevant_dofs,
                 mechanics_communicator);

        for (auto const dof : locally_owned_dofs)
        {
          if (!component_dofs[dof])
            continue;

          auto const point = support_points.find(dof);

          if (point == support_points.end())
            continue;

          unsigned int const displacement_component =
              dof_components[dof];

          if (displacement_component >= dim)
            continue;

          dealii::Tensor<1, dim> r;

          for (unsigned int d = 0; d < dim; ++d)
          {
            r[d] =
                point->second[d] - component_center[d];
          }

          double value = 0.0;

          if (rigid_mode < dim)
          {
            // Translation in coordinate direction rigid_mode.
            if (displacement_component == rigid_mode)
              value = 1.0;
          }
          else if constexpr (dim == 2)
          {
            // One in-plane rotation:
            //
            // omega_z x r = (-r_y, r_x)

            if (displacement_component == 0)
              value = -r[1];
            else if (displacement_component == 1)
              value = r[0];
          }
          else if constexpr (dim == 3)
          {
            unsigned int const rotation_axis =
                rigid_mode - dim;

            // Rotation about x:
            //
            // omega_x x r = (0, -r_z, r_y)
            if (rotation_axis == 0)
            {
              if (displacement_component == 1)
                value = -r[2];
              else if (displacement_component == 2)
                value = r[1];
            }

            // Rotation about y:
            //
            // omega_y x r = (r_z, 0, -r_x)
            else if (rotation_axis == 1)
            {
              if (displacement_component == 0)
                value = r[2];
              else if (displacement_component == 2)
                value = -r[0];
            }

            // Rotation about z:
            //
            // omega_z x r = (-r_y, r_x, 0)
            else if (rotation_axis == 2)
            {
              if (displacement_component == 0)
                value = -r[1];
              else if (displacement_component == 1)
                value = r[0];
            }
          }

          mode[dof] = value;
        }

        mode.compress(dealii::VectorOperation::insert);

        // Modified Gram-Schmidt.
        //
        // Modes from disconnected components have disjoint support.
        // Modes belonging to the same component still need
        // orthogonalization, particularly rotations versus translations.
        for (auto const &existing_mode :
             _floating_rigid_body_modes)
        {
          double const coefficient =
              existing_mode * mode;

          mode.add(-coefficient, existing_mode);
        }

        double const mode_norm = mode.l2_norm();

        if (mode_norm > 1.e-14)
        {
          mode /= mode_norm;

          _floating_rigid_body_modes.push_back(mode);

          ++modes_added_for_component;
        }
      }

      std::cout
          << "Floating mechanical component "
          << floating_component_number
          << ": cells = "
          << component.cell_indices.size()
          << ", rigid-body modes = "
          << modes_added_for_component
          << std::endl;

      ++floating_component_number;
    }

    std::cout
        << "Mechanical active components: "
        << components.size()
        << ", floating components: "
        << floating_component_number
        << ", projected rigid-body modes: "
        << _floating_rigid_body_modes.size()
        << std::endl;
  }
  else
  {
    std::cout
        << "Floating-component nullspace projection is currently "
           "implemented only for one MPI rank with no hanging-node "
           "constraints; projection is disabled for this mechanical "
           "system."
        << std::endl;
  }
}

template <int dim, int n_materials, int p_order, typename MaterialStates,
          typename MemorySpaceType>
dealii::LA::distributed::Vector<double, dealii::MemorySpace::Host>
MechanicalPhysics<dim, n_materials, p_order, MaterialStates,
                  MemorySpaceType>::solve()
{
#ifdef ADAMANTINE_WITH_CALIPER
  CALI_MARK_BEGIN("solve mechanical system");
#endif

  dealii::IndexSet locally_owned_dofs = _dof_handler.locally_owned_dofs();
  dealii::IndexSet locally_relevant_dofs =
      dealii::DoFTools::extract_locally_relevant_dofs(_dof_handler);
#if DEAL_II_VERSION_GTE(9, 7, 0) && defined(DEAL_II_TRILINOS_WITH_TPETRA)
  using TrilinosVectorType = dealii::LinearAlgebra::TpetraWrappers::Vector<
      double, dealii::MemorySpace::Default>;
#else
  using TrilinosVectorType = dealii::TrilinosWrappers::MPI::Vector;
#endif
  TrilinosVectorType displacement(
      locally_owned_dofs, _mechanical_operator->rhs().get_mpi_communicator());
  TrilinosVectorType rhs_device(
      locally_owned_dofs, _mechanical_operator->rhs().get_mpi_communicator());
  dealii::LinearAlgebra::ReadWriteVector<double> rw_vector(locally_owned_dofs);

  rw_vector.import_elements(_mechanical_operator->rhs(),
                            dealii::VectorOperation::insert);
  rhs_device.import_elements(rw_vector, dealii::VectorOperation::insert);
  double const rhs_norm = rhs_device.l2_norm();
  // Solve the mechanical problem assuming that the deformation is elastic
  // TODO check that we are computing only difference of the displacement
  // compared to the previous time step!!
  unsigned int const max_iter = _dof_handler.n_dofs() / 10;
  double const tol = 1e-12 * rhs_norm;
  dealii::SolverControl solver_control(max_iter, tol);
  dealii::SolverCG<TrilinosVectorType> cg(solver_control);

  std::exception_ptr regular_cg_exception;
#ifdef ADAMANTINE_WITH_CALIPER
  CALI_MARK_BEGIN("mechanical regular CG attempt");
#endif
  try
  {
    cg.solve(_mechanical_operator->system_matrix(),
             displacement,
             rhs_device,
             _mechanical_operator->preconditioner());
  }
  catch (dealii::SolverControl::NoConvergence const &)
  {
#ifdef ADAMANTINE_WITH_CALIPER
  CALI_MARK_END("mechanical regular CG attempt");
#endif
    regular_cg_exception = std::current_exception();
  }
  catch (...)
  {
#ifdef ADAMANTINE_WITH_CALIPER
  CALI_MARK_END("mechanical regular CG attempt");
#endif
    throw;
  }
  if (!regular_cg_exception)
  {
#ifdef ADAMANTINE_WITH_CALIPER
  CALI_MARK_END("mechanical regular CG attempt");
#endif
  }
  if (regular_cg_exception)
  {
#ifdef ADAMANTINE_WITH_CALIPER
  CALI_MARK_BEGIN("CG failed - mechanical floating component detection");
#endif
    detect_floating_rigid_body_modes();
#ifdef ADAMANTINE_WITH_CALIPER
  CALI_MARK_END("CG failed - mechanical floating component detection");
#endif
    // If CG failed but there is no floating component, don't pretend the
    // failure was caused by this nullspace. Preserve the original exception.
    if (_floating_rigid_body_modes.empty())
#ifdef ADAMANTINE_WITH_CALIPER
  CALI_MARK_END("solve mechanical system");
#endif
      std::rethrow_exception(regular_cg_exception);
#ifdef ADAMANTINE_WITH_CALIPER
  CALI_MARK_BEGIN("CG failed - mechanical projection fallback setup");
#endif
    std::vector<TrilinosVectorType> floating_rigid_body_modes;
    floating_rigid_body_modes.reserve(_floating_rigid_body_modes.size());
    for (auto const &host_mode : _floating_rigid_body_modes)
    {
      dealii::LinearAlgebra::ReadWriteVector<double> rw_mode(
          locally_owned_dofs);

      rw_mode.import_elements(host_mode,
                              dealii::VectorOperation::insert);

      floating_rigid_body_modes.emplace_back(
          locally_owned_dofs,
          _mechanical_operator->rhs().get_mpi_communicator());

      floating_rigid_body_modes.back().import_elements(
          rw_mode,
          dealii::VectorOperation::insert);
    }
    // The failed ordinary CG has changed displacement.
    // Restart the projected solve from zero.
    displacement = 0.0;
    // rhs_device itself should still contain the original RHS.
    // Project it only now, after failure.
    project_out_nullspace(rhs_device, floating_rigid_body_modes);
    double const projected_rhs_norm = rhs_device.l2_norm();
    double const projected_tol = 1e-12 * projected_rhs_norm;
#ifdef ADAMANTINE_WITH_CALIPER
  CALI_MARK_END("CG failed - mechanical projection fallback setup");
#endif
  using MechanicalOperatorType =
      MechanicalOperator<dim, n_materials, p_order, MaterialStates,
                         MemorySpaceType>;

  using MatrixType =
      typename MechanicalOperatorType::TrilinosMatrixType;

  using PreconditionerType =
      typename MechanicalOperatorType::TrilinosPreconditionerType;

#ifdef ADAMANTINE_WITH_CALIPER
  CALI_MARK_BEGIN("mechanical projected CG retry total");
#endif
  NullspaceProjectedOperator<MatrixType, TrilinosVectorType>
      projected_matrix(
          _mechanical_operator->system_matrix(),
          floating_rigid_body_modes,
          locally_owned_dofs,
          _mechanical_operator->rhs().get_mpi_communicator());
  NullspaceProjectedPreconditioner<PreconditionerType, TrilinosVectorType>
      projected_preconditioner(
          _mechanical_operator->preconditioner(),
          floating_rigid_body_modes,
          locally_owned_dofs,
          _mechanical_operator->rhs().get_mpi_communicator());

  dealii::SolverControl projected_solver_control(max_iter,
                                                 projected_tol);
  dealii::SolverCG<TrilinosVectorType> projected_cg(
      projected_solver_control);
  try
  {
#ifdef ADAMANTINE_WITH_CALIPER
  CALI_MARK_BEGIN("mechanical projected CG solve");
#endif
    projected_cg.solve(projected_matrix,
                       displacement,
                       rhs_device,
                       projected_preconditioner);
#ifdef ADAMANTINE_WITH_CALIPER
  CALI_MARK_END("mechanical projected CG solve");
#endif
    project_out_nullspace(displacement,
                          floating_rigid_body_modes);
  }
  catch (...)
  {
#ifdef ADAMANTINE_WITH_CALIPER
    CALI_MARK_END("mechanical projected CG solve");
    CALI_MARK_END("mechanical projected CG retry total");
#endif
    throw;
  }
#ifdef ADAMANTINE_WITH_CALIPER
  CALI_MARK_END("mechanical projected CG retry total");
#endif
}
  rw_vector.import_elements(displacement, dealii::VectorOperation::insert);
  dealii::LA::distributed::Vector<double, dealii::MemorySpace::Host>
      displacement_host(locally_owned_dofs, locally_relevant_dofs,
                        _mechanical_operator->rhs().get_mpi_communicator());
  displacement_host.import_elements(rw_vector, dealii::VectorOperation::insert);
  _affine_constraints.distribute(displacement_host);

  // Compute the new stress assuming the deformation is elastic.
  // If the stress is under the yield criterion, the deformation is elastic and
  // we are done. Otherwise we need to use the radial return algorithm to
  // compute the plastic deformation.
  dealii::LA::distributed::Vector<double, dealii::MemorySpace::Host>
      incremental_displacement(
          locally_owned_dofs, locally_relevant_dofs,
          _mechanical_operator->rhs().get_mpi_communicator());
  incremental_displacement = displacement_host;
  if (_old_displacement.size() > 0)
  {
    incremental_displacement -= _old_displacement;
  }
  incremental_displacement.update_ghost_values();
  compute_stress(incremental_displacement);

  _old_displacement.swap(displacement_host);

#ifdef ADAMANTINE_WITH_CALIPER
  CALI_MARK_END("solve mechanical system");
#endif

  return _old_displacement;
}

template <int dim, int n_materials, int p_order, typename MaterialStates,
          typename MemorySpaceType>
void MechanicalPhysics<dim, n_materials, p_order, MaterialStates,
                       MemorySpaceType>::
    compute_stress(
        dealii::LA::distributed::Vector<double, dealii::MemorySpace::Host> const
            &displacement)
{
  dealii::hp::FEValues<dim> displacement_hp_fe_values(
      _fe_collection, _q_collection, dealii::update_gradients);
  unsigned int const n_q_points = _q_collection.max_n_quadrature_points();
  std::vector<dealii::SymmetricTensor<2, dim>> strain_tensor(n_q_points);
  std::vector<double> temperature_values(n_q_points);
  const dealii::FEValuesExtractors::Vector displacement_extr(0);
  std::unique_ptr<dealii::hp::FEValues<dim>> temperature_hp_fe_values;
  std::vector<unsigned int> cell_indices;
  // When solving a thermomechanical problem, we create mapping between the
  // active cells of the mechanical problem and the active cells of the thermal
  // problem. Liquid cells are active in the thermal problem but not in the
  // mechanical problem.
  if (!_reference_temperatures.empty())
  {
    _temperature.update_ghost_values();
    temperature_hp_fe_values = std::make_unique<dealii::hp::FEValues<dim>>(
        _thermal_dof_handler->get_fe_collection(), _q_collection,
        dealii::update_values);

    auto &triangulation = _dof_handler.get_triangulation();
    cell_indices.resize(triangulation.n_active_cells());
    unsigned int thermal_cell_index = 0;
    for (auto const &tria_cell :
         triangulation.active_cell_iterators() |
             dealii::IteratorFilters::LocallyOwnedCell())
    {
      dealii::TriaIterator<dealii::DoFCellAccessor<dim, dim, false>>
          temperature_cell(&triangulation, tria_cell->level(),
                           tria_cell->index(), _thermal_dof_handler);
      if (temperature_cell->active_fe_index() == 0)
      {
        dealii::TriaIterator<dealii::DoFCellAccessor<dim, dim, false>>
            displacement_cell(&triangulation, tria_cell->level(),
                              tria_cell->index(), &_dof_handler);
        if (displacement_cell->active_fe_index() == 0)
          cell_indices[displacement_cell->active_cell_index()] =
              thermal_cell_index;
        ++thermal_cell_index;
      }
    }
  }
  unsigned int cell_id = 0;
  for (auto const &cell : _dof_handler.active_cell_iterators())
  {
    if (cell->is_locally_owned() && cell->active_fe_index() == 0)
    {
      // Formulation based on the combined isotropic-kinematic hardening model
      // for J2 plasticity in Chapter 3 of R. Borja, Plasticity: Modeling and
      // Computation, Springer-Verlag, 2013. DOI: 10.1007/978-3-642-38547-6
      //
      // Compute the strain. We get the strain for all the quadrature points at
      // once.
      displacement_hp_fe_values.reinit(cell);
      auto const &fe_values = displacement_hp_fe_values.get_present_fe_values();

      fe_values[displacement_extr].get_function_symmetric_gradients(
          displacement, strain_tensor);

      double reference_temperature = 0.;
      if (!_reference_temperatures.empty())
      {
        auto &triangulation = _dof_handler.get_triangulation();
        dealii::TriaIterator<dealii::DoFCellAccessor<dim, dim, false>>
            temperature_cell(&triangulation, cell->level(), cell->index(),
                             _thermal_dof_handler);
        temperature_hp_fe_values->reinit(temperature_cell);
        auto const &temperature_fe_values =
            temperature_hp_fe_values->get_present_fe_values();
        temperature_fe_values.get_function_values(_temperature,
                                                  temperature_values);
        reference_temperature =
            _has_melted[cell_indices[cell->active_cell_index()]]
                ? _reference_temperatures[temperature_cell->material_id()]
                : _reference_temperatures.back();
      }

      double const lambda = _material_properties.get_mechanical_property(
          cell, StateProperty::lame_first_parameter);
      double const mu = _material_properties.get_mechanical_property(
          cell, StateProperty::lame_second_parameter);
      double const plastic_modulus =
          _material_properties.get_mechanical_property(
              cell, StateProperty::plastic_modulus);
      double const iso_hardening_coef =
          _material_properties.get_mechanical_property(
              cell, StateProperty::isotropic_hardening);
      double const alpha = _material_properties.get_mechanical_property(
          cell, StateProperty::thermal_expansion_coef);
      double const beta = (3. * lambda + 2. * mu) * alpha;
      dealii::SymmetricTensor<4, dim> stiffness_tensor =
          lambda * dealii::outer_product(dealii::unit_symmetric_tensor<dim>(),
                                         dealii::unit_symmetric_tensor<dim>()) +
          2 * mu * dealii::identity_tensor<dim>();
      // Loop over the quadrature points.
      for (auto const q : fe_values.quadrature_point_indices())
      {
        // Compute the trial elastic stress.
        dealii::SymmetricTensor<2, dim> elastic_stress = _stress[cell_id][q];
        elastic_stress += stiffness_tensor * strain_tensor[q];
        if (!_reference_temperatures.empty())
        {
          // Compute the thermal stress due to temperature change since the last
          // time step. The rest of the termal stress in already included in the
          // previous strees.
          double const current_thermal_stress =
              beta * (temperature_values[q] - reference_temperature);
          elastic_stress -=
              (current_thermal_stress - _thermal_stress[cell_id][q]) *
              dealii::unit_symmetric_tensor<dim>();
          _thermal_stress[cell_id][q] = current_thermal_stress;
        }

        auto stress_deviator = dealii::deviator(elastic_stress);
        auto effective_stress = stress_deviator - _back_stress[cell_id][q];
        double const effective_stress_norm = effective_stress.norm();
        if (effective_stress_norm <= _plastic_internal_variable[cell_id][q])
        {
          // The deformation is elastic. We just update the stress with the
          // elastic stress.
          _stress[cell_id][q] = elastic_stress;
        }
        else
        {
          // The deformation is plastic. We need to compute a new stress and
          // update the plastic internal variable and the back stress.
          double plastic_strain_increment =
              (effective_stress_norm - _plastic_internal_variable[cell_id][q]) /
              (2. * mu + plastic_modulus);
          auto plastic_flow_direction =
              effective_stress / effective_stress_norm;
          // Update stress
          _stress[cell_id][q] = elastic_stress - 2. * mu *
                                                     plastic_strain_increment *
                                                     plastic_flow_direction;
          // Update plastic internal variable
          _plastic_internal_variable[cell_id][q] +=
              iso_hardening_coef * plastic_modulus * plastic_strain_increment;
          // Update back stress
          _back_stress[cell_id][q] +=
              (1. - iso_hardening_coef) * plastic_modulus *
              plastic_strain_increment * plastic_flow_direction;
        }
      }
    }
    ++cell_id;
  }
}

} // namespace adamantine

INSTANTIATE_DIM_NMAT_PORDER_MATERIALSTATES_HOST(MechanicalPhysics)
INSTANTIATE_DIM_NMAT_PORDER_MATERIALSTATES_DEVICE(MechanicalPhysics)
