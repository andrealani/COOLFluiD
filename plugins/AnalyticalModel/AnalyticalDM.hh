// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_AnalyticalDM_hh
#define COOLFluiD_Framework_AnalyticalDM_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DomainModel.hh"
#include "Framework/VectorialFunction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace AnalyticalModel {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a model definition
 * described by analytical parametric functions.
 *
 * @author Tiago Quintino
 */
class AnalyticalDM : public Framework::DomainModel {
public: // interface functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  AnalyticalDM(const std::string& name);

  /// Destructor
  virtual ~AnalyticalDM();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// gets the number of topological regions that have definition
  /// @return number of topological regions with definition
  virtual TRidx getNbTopoDefs () const ;

  /// computes the parametric coordinates given the real model coordinates
  /// @pre return parameter must be properly resized
  virtual void computeParamCoord (const TRidx idx, const XVector& coord , PVector& pcoord) const ;

  /// computes the coordinates in the model space given the parametric coordinates
  /// @pre return parameter must be properly resized
  virtual void computeCoord (const TRidx idx, const PVector& pcoord, XVector& coord) const ;

  /// computes the first derivatives in the model space given the parametric coordinates
  /// @pre return parameter must be properly resized
  virtual void compute1stDeriv (const TRidx idx, const PVector& pcoord, std::vector< XVector >& deriv1) const	;

  /// computes the second derivatives in the model space given the parametric coordinates
  /// @pre return parameter must be properly resized
  virtual void compute2ndDeriv (const TRidx idx, const PVector& pcoord, std::vector< XVector >& deriv2) const	;

  /// computes the coordinates, first and second derivatives all in one call
  /// @pre return parameters must be properly resized
	virtual void computeAll (const TRidx idx, const PVector& pcoord, XVector& coord, std::vector< XVector >& deriv1, std::vector< XVector >& deriv2) const ;

  /// Gets the class name
  static std::string getClassName() { return "AnalyticalDM"; }

private: // data

  /// number of topological regions in the model
  Framework::DomainModel::TRidx m_nb_topo;

  /// dimension of the model space
  CFuint m_modeldim;

  /// dimension of the parametic space
  CFuint m_pardim;

  /// a vector of strings for definition of the analytical functions
  /// for the x coordinate
  std::vector<std::string> m_x_func_def;

  /// a vector of strings for definition of the analytical functions
  /// for the y coordinate
  std::vector<std::string> m_y_func_def;

  /// a vector of strings for definition of the analytical functions
  /// for the z coordinate
  std::vector<std::string> m_z_func_def;

  /// the definition of the coordinate functions
  std::vector< Framework::VectorialFunction > m_func;

  /// a vector of strings for definition of the 1st derivatives analytical functions
  /// for the x coordinate
  std::vector<std::string> m_x_dfunc_def;

  /// a vector of strings for definition of the 1st derivatives analytical functions
  /// for the y coordinate
  std::vector<std::string> m_y_dfunc_def;

  /// a vector of strings for definition of the 1st derivatives analytical functions
  /// for the z coordinate
  std::vector<std::string> m_z_dfunc_def;

  /// the definition of the 1st derivative functions
  std::vector< Framework::VectorialFunction > m_dfunc;

}; // end of class AnalyticalDM

//////////////////////////////////////////////////////////////////////////////

  } // namespace AnalyticalModel

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_AnalyticalDM_hh
