#ifndef COOLFluiD_Numerics_SpectralFV_NSFaceDiffusiveFluxLocalApproach_hh
#define COOLFluiD_Numerics_SpectralFV_NSFaceDiffusiveFluxLocalApproach_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "NavierStokes/NavierStokesVarSet.hh"

#include "SpectralFV/FaceDiffusiveFluxLocalApproach.hh"
#include "SpectralFV/SpectralFVMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the averaged diffusive flux on a face,
 * based on the `local DG' approach, for the Navier-Stokes equations
 *
 * @author Kris Van den Abeele
 */
class NSFaceDiffusiveFluxLocalApproach : public FaceDiffusiveFluxLocalApproach {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFVMethodData,NSFaceDiffusiveFluxLocalApproach > PROVIDER;

public:  // methods

  /// Constructor
  NSFaceDiffusiveFluxLocalApproach(const std::string& name);

  /// Destructor
  ~NSFaceDiffusiveFluxLocalApproach();

  /// Compute averaged gradient variables in a series of flux points, from left and right states and face normal
  virtual std::vector< RealVector >& computeAvgGradVars(std::vector< RealVector* >& lStates,
                                                        std::vector< RealVector* >& rStates,
                                                        const CFuint nbrFlxPnts);

  /// Compute the averaged diffusive flux in a series of flux points,
  /// from left and right gradients and states and a normal vector
  virtual std::vector< RealVector >& computeDiffFlux(std::vector< std::vector< RealVector* >* >& lGrads,
                                                     std::vector< std::vector< RealVector* >* >& rGrads,
                                                     std::vector< RealVector* >& lStates,
                                                     std::vector< RealVector* >& rStates,
                                                     const CFreal& lCellVol,
                                                     const CFreal& rCellVol,
                                                     const CFreal& faceSurf,
                                                     const std::vector< RealVector >& normal,
                                                     const CFuint nbrFlxPnts);

  /// Gets the Class name
  static std::string getClassName()
  {
    return "NSFaceDiffusiveFluxLocalApproach";
  }

  /// Setup private data
  virtual void setup();

  /// Unsetup private data
  virtual void unsetup();

protected: // data

  /// Navier-Stokes variable set
  Common::SafePtr< Physics::NavierStokes::NavierStokesVarSet > m_navierStokesVarSet;

  /// gradient variables in the flux points in left cell
  RealMatrix m_flxPntLGradVars;

  /// gradient variables in the flux points in right cell
  RealMatrix m_flxPntRGradVars;

}; // class NSFaceDiffusiveFluxLocalApproach

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFV_NSFaceDiffusiveFluxLocalApproach_hh

