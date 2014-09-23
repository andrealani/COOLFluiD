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

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a model definition
/// described by analytical parametric functions.
/// @author Tiago Quintino
class Framework_API NullDomainModel : public Framework::DomainModel {
public: // interface functions

  /// Constructor
  NullDomainModel(const std::string& name);

  /// Destructor
  virtual ~NullDomainModel();

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
  virtual void compute1stDeriv (const TRidx idx, const PVector& pcoord, std::vector< XVector >& deriv1) const  ;

  /// computes the second derivatives in the model space given the parametric coordinates
  /// @pre return parameter must be properly resized
  virtual void compute2ndDeriv (const TRidx idx, const PVector& pcoord, std::vector< XVector >& deriv2) const  ;

  /// computes the coordinates, first and second derivatives all in one call
  /// @pre return parameters must be properly resized
virtual void computeAll (const TRidx idx, const PVector& pcoord, XVector& coord, std::vector< XVector >& deriv1, std::vector< XVector >& deriv2) const ;

  /// Gets the class name
  static std::string getClassName() { return "NullDomainModel"; }

}; // end of class NullDomainModel

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_AnalyticalDM_hh
