/* Copyright (c) 2020 - 2024, the adamantine authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#include <DonutHeatSource.hh>
#include <instantiation.hh>
#include <types.hh>
#include <cmath> // For mathematical functions

namespace adamantine
{

// Constructor: Initializes the donut-shaped heat source for AFSD (Ti64)
template <int dim>
DonutHeatSource<dim>::DonutHeatSource(boost::property_tree::ptree const &database)
    : HeatSource<dim>(database),
      _current_angle(0.0)
{
  // Retrieve the inner radius and diameter of the donut shape from the input file
  double inner_radius = database.get<double>("inner_radius");
  _inner_radius_squared = std::pow(inner_radius, 2);
  _diameter = database.get<double>("diameter");
  // Retrieve rotational parameters for AFSD
  _rotation_speed = database.get<double>("rotation_speed", 0.0); // Default to 0 if not provided
  _advancing_heat_modifier = database.get<double>("advancing_heat_modifier", 1.0);
  _retreating_heat_modifier = database.get<double>("retreating_heat_modifier", 1.0);
  // Retrieve the material thickness from the input file
  _material_thickness = database.get<double>("depth");
  // Define a constant factor k for the heat input calculation
  _efficiency = database.get<double>("absorption_efficiency", 0.8); // Default to 0.8, adjust based on calibration
  _force = database.get<double>("force", 4900.0); // Default to 1.0, adjust based on calibration
  _sigma = database.get<double>("heat_distribution_sigma", _diameter*4.0); // Control heat spread
}

// Update the current time step and calculate heat input based on simplified factors
template <int dim>
void DonutHeatSource<dim>::update_time(double time)
{
  // Update the position of the beam (tool center) based on the scan path at the current time
  _beam_center = this->_scan_path.value(time);

  // Calculate the current angle based on rotation speed and time for rotating tool behavior
  _current_angle = std::fmod(_rotation_speed * time, 2.0 * dealii::numbers::PI);

  // Calculate tool radius (half the diameter)
  double tool_radius = _diameter / 2.0;

  // Calculate angular velocity (in rad/s)
  double angular_velocity = _rotation_speed * 2.0 * dealii::numbers::PI / 60.0; // Convert RPM to rad/s

  // Calculate torque (Force * Radius)
  double torque = _force * tool_radius;

  // Calculate power (Torque * Angular Velocity)
  double power = torque * angular_velocity;

  // Calculate heat input based on power and efficiency
  _alpha = _efficiency * power;
  std::cout << "Heat input value: " << _alpha << std::endl;
}

// This method calculates the heat source value at a given point in space and time.
template <int dim>
double DonutHeatSource<dim>::value(dealii::Point<dim> const &point, double const height) const
{
  // z-axis adjustment based on the height of the material
  double const z = point[axis<dim>::z] - height;

  // If the point is below the tool depth, no heat contribution
  if ((z + this->_beam.depth) < 0.)
  {
    return 0.;
  }
  else
  {
    // Calculate the radial distance from the beam center in the x-y plane
    double xpy_squared = std::pow(point[axis<dim>::x] - _beam_center[axis<dim>::x], 2);
    if (dim == 3)
    {
      xpy_squared += std::pow(point[axis<dim>::y] - _beam_center[axis<dim>::y], 2);
    }

    // Exclude points within the inner radius of the donut-shaped heat source (center of tool)
    if (xpy_squared < _inner_radius_squared)
    {
      return 0.;
    }

    // Calculate the angle of the point relative to the current rotation angle of the tool
    double point_angle = std::atan2(point[axis<dim>::y], point[axis<dim>::x]);
    double relative_angle = std::fmod(point_angle - _current_angle, 2.0 * dealii::numbers::PI);

    // Determine if the point is on the advancing side or retreating side
    double heat_modifier = (relative_angle >= -dealii::numbers::PI / 2.0 && relative_angle <= dealii::numbers::PI / 2.0)
                               ? _advancing_heat_modifier
                               : _retreating_heat_modifier;
    std::cout << "r " << xpy_squared << std::endl;
    std::cout << "denom " << 2.0 * std::pow(_sigma, 2) << std::endl;
    std::cout << "exp " << std::exp(-xpy_squared / (2.0 * std::pow(_sigma, 2))) << std::endl;
    // Apply a Gaussian radial decay for FSW heat distribution
    double heat_source =1000* _alpha * heat_modifier * std::exp(-xpy_squared / (2.0 * std::pow(_sigma, 2)));
    std::cout << "Heat source value: " << heat_source << std::endl;
    return heat_source;
  }
}

} // namespace adamantine

INSTANTIATE_DIM(DonutHeatSource)
