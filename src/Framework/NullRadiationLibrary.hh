// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_Framework_NullRadiationLibrary_hh
#define COOLFluiD_Framework_NullRadiationLibrary_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/RadiationLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
        
//////////////////////////////////////////////////////////////////////////////

/// This class represents the interface for a catalycity library
/// @author Andrea Lani
class Framework_API NullRadiationLibrary : public Framework::RadiationLibrary {
public:
  
  /// Constructor without arguments
  NullRadiationLibrary(const std::string& name);
  
  /// Default destructor
  virtual ~NullRadiationLibrary();
  
  /// set up private data
  virtual void setup();
  
  /// Run the radiation code on the stagnation line
  virtual void runOnStagnationLine(Common::SafePtr<std::vector<CFuint> > stagnationLineCells,
				   Framework::ProxyDofIterator<CFreal>* pstates,
				   CFreal* qrad);
  
  /// Run the radiation code on a structured mesh
  virtual void runOnStructuredMesh(const std::vector<std::vector<CFuint>* >& meshByLine,
				   Framework::ProxyDofIterator<CFreal>* pstates,
				   CFreal* qrad);
  
  /// Compute radiative properties
  virtual void computeProperties(Framework::ProxyDofIterator<CFreal>* pstates,
				 RealMatrix& data, CFuint iWavRange);
  
}; // end of class NullRadiationLibrary
    
//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NullRadiationLibrary_hh
