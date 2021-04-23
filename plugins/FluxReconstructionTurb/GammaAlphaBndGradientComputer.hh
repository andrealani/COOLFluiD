// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_GammaAlphaBndGradientComputer_hh
#define COOLFluiD_FluxReconstructionMethod_GammaAlphaBndGradientComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionNavierStokes/NSBndGradientComputer.hh"
#include "NavierStokes/NavierStokesVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Daughterclass of NSBndGradientComputer, needed to 
 * calculate the bnd gradients for implicit schemes for gamma-alpha
 * 
 * @author Ray Vandenhoeck
 */
class GammaAlphaBndGradientComputer : public NSBndGradientComputer {

public: // functions

  /// Constructor
  explicit GammaAlphaBndGradientComputer(const std::string& name);

  /// Destructor
  virtual ~GammaAlphaBndGradientComputer() {}
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();

protected:
    
  /// compute the states and ghost states in the flx pnts
  virtual void computeFlxPntStates();
  
protected: //data
    
  /// storage of the volumes
  Framework::DataSocketSink<CFreal> socket_volumes;
  
  /// storage for the normals in the solution points
  Framework::DataSocketSink< CFreal > socket_solPntNormals;
    
  /// dummy vector for gradients
  std::vector< RealVector* > m_gradDummy;
  
  /// u gradient
  RealVector m_uGrad;
  
  /// v gradient
  RealVector m_vGrad;
  
  /// w gradient
  RealVector m_wGrad;
  
  /// flux projection vectors in solution points for disc flux
  std::vector< std::vector< RealVector > > m_cellFluxProjVects;
  
  /// dependencies of solution pnts on sol pnts
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_solSolDep;
  
  /// nbr of sol pnts on which a flx pnt is dependent
  CFuint m_nbrSolSolDep;
  
  /// coefs to compute the derivative of the states in the sol pnts
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > m_solPolyDerivAtSolPnts;
  
  /// matrix to store the state terms needed for the gradients (p, u, v, T) inside element
  RealMatrix m_tempGradTermIntCell;
  
  /// element states within an element in the correct format
  std::vector< RealVector* > m_tempStatesIntCell;
    
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_GammaAlphaBndGradientComputer_hh

