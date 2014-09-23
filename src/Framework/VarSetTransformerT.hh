// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_VarSetTransformerT_hh
#define COOLFluiD_Framework_VarSetTransformerT_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/COOLFluiD.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the basic interface for a variable transformer.
/// It accepts 3 template parameters, corresponding to 3 variable sets
/// The last parameter should be NOTYPE by default: when specified, VarSetTransformerT 
/// indicates a matrix transformer, otherwise it behaves as a vector transformer.
/// @author Andrea Lani

template <typename FROM, typename TO, typename IN>     
class VarSetTransformerT {};
    
//////////////////////////////////////////////////////////////////////////////

/// This class implements a partial specialization to force users to implement
/// transformations with FROM different from TO. This serves as default 
/// implementation for cases when FROM is equal to TO (whatever IN).
/// @author Andrea Lani

template <typename FROM, typename IN>     
class VarSetTransformerT<FROM,FROM,IN> {
public:
  
  /// Constructor
  HOST_DEVICE VarSetTransformerT(typename FROM::PTERM::template DeviceConfigOptions<NOTYPE>* dco) {}
  
  /// Default constructor
  HOST_DEVICE VarSetTransformerT() {}
  
  /// Destructor
  HOST_DEVICE virtual ~VarSetTransformerT() {}
  
  /// set the model data
  HOST_DEVICE void setModelData(typename FROM::PTERM::template DeviceConfigOptions<NOTYPE>* dco) {}
  
  /// Transform a state into another one (no transformation by default)
  HOST_DEVICE void transform(const CFreal *const state, CFreal *const result)
  {
    for (CFuint i = 0; i < FROM::NBEQS; ++i) {result[i] = state[i];}
  }
};
    
//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_VarSetTransformerT_hh
