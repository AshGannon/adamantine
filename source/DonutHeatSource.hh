/* Copyright (c) 2020 - 2024, the adamantine authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */
#ifndef DONUT_HEAT_SOURCE_HH
#define DONUT_HEAT_SOURCE_HH
#include <HeatSource.hh>

#include <limits>

namespace adamantine
{
 /**
 * A derived class from HeatSource for a donut shape to be used for
 * friction stir welding.
*/
template <int dim>
class DonutHeatSource final : public HeatSource<dim>
{
public:
  /**
   * Constructor.
   * \param[in] database requires the following entries:
   *   - <B>absorption_efficiency</B>: double in \f$[0,1]\f$
   *   - <B>depth</B>: double in \f$[0,\infty)\f$
   *   - <B>diameter</B>: double in \f$[0,\infty)\f$
   *   - <B>inner_radius</B>: double in \f$[0,\infty)\f$
   *   - <B>max_power</B>: double in \f$[0, \infty)\f$
   *   - <B>advancing_heat_modifier</B>: optional, double with default value \f$1.0\f$
   *   - <B>retreating_heat_modifier</B>: optional, double with default value \f$1.0\f$
   *   - <B>rotation_speed</B>: optional, angular velocity in radians per second, default is \f$0.0\f$
   *   - <B>input_file</B>: name of the file that contains the scan path
   *     segments
   */
  DonutHeatSource(boost::property_tree::ptree const &database);

  /**
   * Set the time variable.
   */
  void update_time(double time) final;

  /**
   * Returns the value of a donut-shaped heat source at a specified point and
   * time.
   */
  double value(dealii::Point<dim> const &point,
               double const height) const final;
  double compute_flow_stress(double temperature) const;

private:
  dealii::Point<3> _beam_center;
  double _alpha = std::numeric_limits<double>::signaling_NaN();
  double const _pi_over_3_to_1p5 = std::pow(dealii::numbers::PI / 3.0, 1.5);
  double _inner_radius_squared;
  double _diameter;
  // Member variables for rotational dynamics of the heat source
  double _rotation_speed; // Angular velocity (radians per second)
  double _current_angle;  // Current angle based on time
  //account for the asymmetry in heat generation between the advancing and retreating sides of the tool.
  double _advancing_heat_modifier;
  double _retreating_heat_modifier;
  double _material_thickness;
  double _efficiency; // Constant factor for heat input calculation
  double _force;
  double _sigma;
};
} // namespace adamantine

#endif



