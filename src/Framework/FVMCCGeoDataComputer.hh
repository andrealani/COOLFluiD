// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_FVMCCGeoDataComputer_hh
#define COOLFluiD_Framework_FVMCCGeoDataComputer_hh

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

/// This class offers a basic interface for (re-)computing geometric data 
/// (normals, volumes, etc.) for a cell-centered Finite Volume discretization, 
/// assuming that all data arrays have been already resized correctly
/// @author Andrea Lani
template < typename METHODDATA >
class FVMCCGeoDataComputer : public Framework::GeoDataComputer<METHODDATA> {
public:

  /// Constructor
  FVMCCGeoDataComputer(const std::string& name);

  /// Default destructor
  virtual ~FVMCCGeoDataComputer();
  
  /// Set private data that will be used during the computation
  virtual void setup();
  
  /// Compute the geometric data and fill pre-allocated data arrays 
  virtual void compute();
  
  /// Compute the nodes at intermediate time for the cell centers
  virtual void modifyOffMeshNodes();
  
  /// Unsetup private data
  virtual void unsetup();
   
  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
protected: // helper functions
  
  /// Compute the normals data
  virtual void computeNormalsData();

  /// Compute the cell volumes
  virtual void computeCellVolumes();
  
  /// Compute the face centers
  virtual void computeFaceCenters();
  
protected:
  
  /// storage of face normals
  Framework::DataSocketSink<CFreal> socket_normals;
  
  /// storage of face areas
  Framework::DataSocketSink<CFreal> socket_faceAreas;
  
  /// storage of face centroids
  Framework::DataSocketSink<CFreal> socket_faceCenters;
  
  /// IDs corresponding to the cell for which the normal point outward
  Framework::DataSocketSink<CFint> socket_isOutward;
  
  /// storage of ghost states
  /// this state list has ownership on the states
  Framework::DataSocketSink<Framework::State*> socket_gstates;
  
  /// storage of the cell volumes
  Framework::DataSocketSink<CFreal> socket_volumes;
  
}; // end of class FVMCCGeoDataComputer
      
//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "Framework/FVMCCGeoDataComputer.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_FVMCCGeoDataComputer_hh
