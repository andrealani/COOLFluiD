#ifndef COOLFluiD_Numerics_SpectralFD_NSFaceDiffusiveFluxIPApproach_hh
#define COOLFluiD_Numerics_SpectralFD_NSFaceDiffusiveFluxIPApproach_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "NavierStokes/NavierStokesVarSet.hh"

#include "SpectralFD/BCStateComputer.hh"
#include "SpectralFD/FaceDiffusiveFluxIPApproach.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the averaged diffusive flux on a face,
 * based on the `Interior Penalty' approach, for the Navier-Stokes equations
 *
 * @author Kris Van den Abeele
 */
class NSFaceDiffusiveFluxIPApproach : public FaceDiffusiveFluxIPApproach {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,NSFaceDiffusiveFluxIPApproach > PROVIDER;

public:  // methods

  /// Constructor
  NSFaceDiffusiveFluxIPApproach(const std::string& name);

  /// Destructor
  ~NSFaceDiffusiveFluxIPApproach();

  /// Compute the averaged diffusive flux in a series of flux points,
  /// from left and right gradients and states and a normal vector
  virtual std::vector< RealVector >& computeDiffFlux(std::vector< std::vector< RealVector* >* >& lGrads,
                                                     std::vector< std::vector< RealVector* >* >& rGrads,
                                                     std::vector< RealVector* >& lStates,
                                                     std::vector< RealVector* >& rStates,
                                                     const std::vector< CFreal >& faceInvCharLength,
                                                     const std::vector< RealVector >& normal,
                                                     const CFuint nbrFlxPnts);

  /// Compute the diffusive flux in a series of flux points on a boundary face
//   std::vector< RealVector >& computeBndDiffFlux(const std::vector< std::vector< RealVector* > >& intGrads,
//                                                 const std::vector< RealVector* >& intStates,
//                                                 const std::vector< RealVector* >& ghostStates,
//                                                 const std::vector< CFreal >& faceInvCharLength,
//                                                 const std::vector< RealVector >& normal,
//                                                 const std::vector< RealVector >& flxPntCoords,
//                                                 const CFuint nbrFlxPnts);

  /// Set the boundary condition state computer
//   void setBCStateComputer(Common::SafePtr< BCStateComputer > bcStateComputer)
//   {
//     m_bcStateComputer = bcStateComputer;
//   }

  /// Gets the Class name
  static std::string getClassName()
  {
    return "NSFaceDiffusiveFluxIPApproach";
  }

  /// Setup private data
  virtual void setup();

  /// Unsetup private data
  virtual void unsetup();

protected: // data

  /// Navier-Stokes variable set
  Common::SafePtr< Physics::NavierStokes::NavierStokesVarSet > m_navierStokesVarSet;

  /// Strategy that sets the ghost states corresponding to the boundary condition
//   Common::SafePtr< BCStateComputer> m_bcStateComputer;

  /// gradient variables in the flux points in left cell
  RealMatrix m_flxPntLGradVars;

  /// gradient variables in the flux points in right cell
  RealMatrix m_flxPntRGradVars;

  /// states in flux points
  std::vector<RealVector*> m_avgStates;

  /// gradient variable gradients in flux points
  std::vector< std::vector<RealVector*> > m_stateGradients;

  /// gradient variable gradients in flux points
  std::vector< std::vector<RealVector*> > m_gradVarGradients;

  /// ghost gradient variable gradients in flux points
//   std::vector< std::vector<RealVector*> > m_ghostGradVarGradients;

}; // class NSFaceDiffusiveFluxIPApproach

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_NSFaceDiffusiveFluxIPApproach_hh

