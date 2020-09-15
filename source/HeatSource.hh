/* Copyright (c) 2016 - 2020, the adamantine authors.
 *
 * This file is subject to the Modified BSD License and may not be distributed
 * without copyright and license information. Please refer to the file LICENSE
 * for the text and further information on this license.
 */

#ifndef HEAT_SOURCE_HH
#define HEAT_SOURCE_HH

#include <utils.hh>

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>

#include <iostream>
#include <istream>
#include <vector>

namespace adamantine
{
/**
 * This enum distinguishes between the two types of scan path segments.
 */
enum class ScanPathSegmentType
{
  line,
  point
};

/**
 * This structure stores the relevant information for a single segment. The scan
 * path input file distingishes between spots and lines, but when we know the
 * end time and end location, spots become lines with a start point equal to its
 * end point. Everything one needs can be determined from these three quantities
 * (and the segment info from the preceding segment) but in the future it might
 * be worth adding in some redundant information like start time/point and
 * velocity.
 */
struct ScanPathSegment
{
  double end_time;            // Unit: seconds
  double power_modifier;      // Dimensionless
  dealii::Point<2> end_point; // Unit: m (NOTE: converted from mm in the file)
};

/**
 * This structure stores all the physical properties necessary to define an
 * heat source.
 */
struct HeatSourceProperties
{
public:
  /**
   * Absolute penetration of the electron beam into the material where 99% of
   * the beam energy is absorbed.
   */
  double depth;
  /**
   * Energy conversion efficiency on the surface.
   */
  double absorption_efficiency;

  /**
   * Square of the beam radius.
   */
  double radius_squared;
  /**
   * Maximum power of the beam.
   */
  double max_power;
};

namespace internal
{
class BeamCenter
{
public:
  BeamCenter();

  void load_segment_list(std::vector<ScanPathSegment> segment_list);

  double value(dealii::Point<1> const &time,
               unsigned int const component = 0) const;

  void rewind_time();

  void save_time();

private:
  mutable unsigned int _current_segment;
  unsigned int _saved_segment;
  mutable dealii::Point<1> _current_time;
  dealii::Point<1> _saved_time;
  std::vector<ScanPathSegment> _segment_list;
};
} // namespace internal

/**
 * Forward declaration of the tester friend class to HeatSource.
 */
class HeatSourceTester;

/**
 * This class describes the evolution of an electron beam source.
 */
template <int dim>
class HeatSource : public dealii::Function<dim>
{
  friend class HeatSourceTester;

public:
  /**
   * Constructor.
   * \param[in] database requires the following entries:
   *   - <B>energy_conversion_efficiency</B>: double in \f$[0,1]\f$
   *   - <B>control_efficiency</b>: double in \f$[0,1]\f$
   *   - <B>depth</B>: double in \f$[0,\infty)\f$
   *   - <B>diameter</B>: double in \f$[0,\infty)\f$
   *   - <B>max_power</B>: double in \f$[0, \infty)\f$ [optional: if not
   *   defined, <i>current</i> and <i>voltage</i> need to be defined]
   *   - <B>current</B>: double in \f$[0, \infty)\f$ [optional: if defined
   *   <i>voltage</i> should be defined too, if not defined <i>max_power</i>
   *   should be defined]
   *   - <B>voltage</B>: double in \f$[0,\infty)\f$ [optional: if defined
   *   <i>current</i> should be defined too, if not defined <i>max_power</i>
   *   should be defined]
   *   - <B>input_file</B>: name of the csv file that contains the successive
   *   position of the electron beam [optional: if not defined then
   *   <i>abscissa</i> and, in 3D, <i>ordinate</i> need to be defined]
   *   - <B>delimiter</B>: delimiting character used in <i>input_file</i>
   *   [required if <i>input_file</i> is defined]
   *   - <B>abscissa</B>: string, abscissa of the beam as a function of time
   *   (e.g. "(t-1) * (t-2)") [optional: need to be defined if <i>input_file</i>
   *   is not defined]
   *   - <B>ordinate</B>: string, ordinate of the beam as a function of time
   *   [required only for three dimensional calculation and if <i>input_file</i>
   *   is not defined]
   */
  HeatSource(boost::property_tree::ptree const &database);

  /**
   * Set the maximum height of the domain. This is the height at which the
   * electron beam penetrate the material.
   */
  void set_max_height(double height);

  /**
   * Compute the heat source at a given point at the current time.
   */
  double value(dealii::Point<dim> const &point,
               unsigned int const component = 0) const override;

  /**
   * Reset the current time and the position to the last saved state.
   */
  void rewind_time();

  /**
   * Save the current time and the position in the list of successive positions
   * of the beam.
   */
  void save_time();

private:
  /**
   * Height of the domain.
   */
  double _max_height;
  /**
   * Structure of the physical properties of the heat source.
   */
  HeatSourceProperties _beam;

  /**
   * The list of segments in the scan path.
   */
  std::vector<ScanPathSegment> _segment_list;

  /**
   * This function reads the scan path file and populates the vector of
   * ScanPathSegments.
   */
  void parse_scan_path(std::string scan_path_file);

  internal::BeamCenter _beam_center;

};

template <int dim>
inline void HeatSource<dim>::set_max_height(double height)
{
  _max_height = height;
}
} // namespace adamantine

#endif
