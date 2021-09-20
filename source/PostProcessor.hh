/* Copyright (c) 2016 - 2021, the adamantine authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#ifndef POST_PROCESSOR_HH
#define POST_PROCESSOR_HH

#include "MaterialProperty.hh"
#include "types.hh"

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/numerics/data_out.h>

#include <boost/property_tree/ptree.hpp>

namespace adamantine
{
/**
 * This class outputs the results using the vtu format.
 */
template <int dim>
class PostProcessor
{
public:
  /**
   * Constructor for non-ensemble simulations.
   * \param database requires the following entries:
   *   - <B>database</B>: boost::property_tree::ptree
   *   - <B>dof_handler</B>: dealii::DoFHandler<dim>
   * material_properties
   */
  PostProcessor(MPI_Comm const &communicator,
                boost::property_tree::ptree const &database,
                dealii::DoFHandler<dim> &dof_handler);
  /**
   * Constructor for ensemble simulations.
   * \param database requires the following entries:
   *   - <B>database</B>: boost::property_tree::ptree
   *   - <B>dof_handler</B>: dealii::DoFHandler<dim>
   *   - <B>ensemble_member_index</B>: int
   */
  PostProcessor(MPI_Comm const &communicator,
                boost::property_tree::ptree const &database,
                dealii::DoFHandler<dim> &dof_handler,
                int ensemble_member_index);

  /**
   * Output the different vtu and pvtu files.
   */
  void output_pvtu(
      unsigned int cycle, unsigned int n_time_step, double time,
      dealii::LA::distributed::Vector<double> const &solution,
      std::array<dealii::LA::distributed::Vector<double>,
                 static_cast<unsigned int>(MaterialState::SIZE)> const &state,
      dealii::DoFHandler<dim> const &material_dof_handler);

  /**
   * Output the pvd file for Paraview.
   */
  void output_pvd();

private:
  /**
   * MPI communicator.
   */
  MPI_Comm _communicator;
  /**
   * Prefix of the different output files.
   */
  std::string _filename_prefix;
  /**
   * Vector of pair of time and pvtu file.
   */
  std::vector<std::pair<double, std::string>> _times_filenames;
  /**
   * DataOut associated with the post-processing.
   */
  dealii::DataOut<dim> _data_out;
  /**
   * DoFHandler associated with the simulation.
   */
  dealii::DoFHandler<dim> &_dof_handler;
};
} // namespace adamantine
#endif
