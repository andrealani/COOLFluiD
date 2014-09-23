// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullGeoDataComputer_hh
#define COOLFluiD_Framework_NullGeoDataComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GeoDataComputer.hh"
#include "Framework/BaseDataSocketSink.hh"
#include "Framework/MethodStrategy.hh"
#include "Framework/Node.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    
//////////////////////////////////////////////////////////////////////////////

/// This class offers a basic interface for (re-)computing algorithm-dependent
/// geometric data (normals, volumes, etc.), assuming that all data arrays have
/// been already resized correctly
/// @author Andrea Lani
template < typename METHODDATA >
class NullGeoDataComputer : public Framework::GeoDataComputer<METHODDATA> {
public:

  /// Constructor
  NullGeoDataComputer(const std::string& name);

  /// Default destructor
  virtual ~NullGeoDataComputer();
  
  /// Set private data that will be used during the computation
  virtual void setup();
  
  /// Compute the geometric data and fill pre-allocated data arrays 
  virtual void compute();
   
  /// Unsetup private data
  virtual void unsetup();
  
  /// Compute the nodes at intermediate time for the cell centers
  virtual void modifyOffMeshNodes();
  
}; // end of class NullGeoDataComputer
      
//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "Framework/NullGeoDataComputer.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullGeoDataComputer_hh
