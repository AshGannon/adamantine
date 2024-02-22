/* Copyright (c) 2016 - 2024, the adamantine authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#include <ThermalPhysics.templates.hh>
#include <instantiation.hh>

INSTANTIATE_DIM_FEDEGREE_QUAD_HOST(TUPLE(ThermalPhysics))

namespace adamantine
{
template class ThermalPhysics<2, 1, dealii::MemorySpace::Default,
                              dealii::QGauss<1>>;
template class ThermalPhysics<2, 2, dealii::MemorySpace::Default,
                              dealii::QGauss<1>>;
template class ThermalPhysics<2, 3, dealii::MemorySpace::Default,
                              dealii::QGauss<1>>;
template class ThermalPhysics<2, 4, dealii::MemorySpace::Default,
                              dealii::QGauss<1>>;
template class ThermalPhysics<2, 5, dealii::MemorySpace::Default,
                              dealii::QGauss<1>>;

template class ThermalPhysics<3, 1, dealii::MemorySpace::Default,
                              dealii::QGauss<1>>;
template class ThermalPhysics<3, 2, dealii::MemorySpace::Default,
                              dealii::QGauss<1>>;
template class ThermalPhysics<3, 3, dealii::MemorySpace::Default,
                              dealii::QGauss<1>>;
template class ThermalPhysics<3, 4, dealii::MemorySpace::Default,
                              dealii::QGauss<1>>;
template class ThermalPhysics<3, 5, dealii::MemorySpace::Default,
                              dealii::QGauss<1>>;

template class ThermalPhysics<2, 1, dealii::MemorySpace::Default,
                              dealii::QGaussLobatto<1>>;
template class ThermalPhysics<2, 2, dealii::MemorySpace::Default,
                              dealii::QGaussLobatto<1>>;
template class ThermalPhysics<2, 3, dealii::MemorySpace::Default,
                              dealii::QGaussLobatto<1>>;
template class ThermalPhysics<2, 4, dealii::MemorySpace::Default,
                              dealii::QGaussLobatto<1>>;
template class ThermalPhysics<2, 5, dealii::MemorySpace::Default,
                              dealii::QGaussLobatto<1>>;

template class ThermalPhysics<3, 1, dealii::MemorySpace::Default,
                              dealii::QGaussLobatto<1>>;
template class ThermalPhysics<3, 2, dealii::MemorySpace::Default,
                              dealii::QGaussLobatto<1>>;
template class ThermalPhysics<3, 3, dealii::MemorySpace::Default,
                              dealii::QGaussLobatto<1>>;
template class ThermalPhysics<3, 4, dealii::MemorySpace::Default,
                              dealii::QGaussLobatto<1>>;
template class ThermalPhysics<3, 5, dealii::MemorySpace::Default,
                              dealii::QGaussLobatto<1>>;
} // namespace adamantine
