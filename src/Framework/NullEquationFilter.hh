// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullEquationFilter_hh
#define COOLFluiD_Framework_NullEquationFilter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/EquationFilter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////
    
/// This class represents the a Null equation filter
/// @author Andrea Lani
template < typename METHODDATA >
class NullEquationFilter : public Framework::EquationFilter<METHODDATA> {
public:

  /// Constructor
  NullEquationFilter(const std::string& name) : 
    Framework::EquationFilter<METHODDATA>(name)
  {    
  }

  /// Default destructor
  ~NullEquationFilter()
  {
  }
  
  /// Set up private data to prepare the simulation
  virtual void setup()
  {  
    Framework::EquationFilter<METHODDATA>::setup();
    CFLog(VERBOSE,  "NullEquationFilter::setup() called!\n");
  }
  
  /// Unsetup up private data to prepare the simulation
  virtual void unsetup()
  { 
    Framework::EquationFilter<METHODDATA>::unsetup();
    CFLog(VERBOSE,  "NullEquationFilter::unsetup() called!\n");
  }
  
  /// Reset all data before starting a new iteration
  virtual void reset() 
  { 
    CFLog(DEBUG_MED,  "NullEquationFilter::reset() called!\n");
  }
  
  /// Filter the equation on the current geometric entity
  virtual bool filterOnGeo(GeometricEntity *const geo) 
  {
    CFLog(DEBUG_MED,  "NullEquationFilter::filterOnGeo() called!\n");
    return true;
  }
  
}; // end of class NullEquationFilter

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullEquationFilter_hh
