#ifndef COOLFluiD_FluxReconstructionEntropyMHD_BCSuperOutletEntropyMHD_hh
#define COOLFluiD_FluxReconstructionEntropyMHD_BCSuperOutletEntropyMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * Supersonic outlet BC for 3D entropy-variable MHD (EntropyMHD).
 *
 * Zero-gradient on all state variables except phi (absorbing sign flip).
 *
 * @author Rayan Dhib
 */
class BCSuperOutletEntropyMHD : public BCStateComputer {
public:

  BCSuperOutletEntropyMHD(const std::string& name);
  virtual ~BCSuperOutletEntropyMHD();

  virtual void setup();

  void computeGhostStates(const std::vector<Framework::State*>& intStates,
                           std::vector<Framework::State*>& ghostStates,
                           const std::vector<RealVector>& normals,
                           const std::vector<RealVector>& coords);

  void computeGhostGradients(const std::vector<std::vector<RealVector*> >& intGrads,
                              std::vector<std::vector<RealVector*> >& ghostGrads,
                              const std::vector<RealVector>& normals,
                              const std::vector<RealVector>& coords);

}; // end of class BCSuperOutletEntropyMHD

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionEntropyMHD_BCSuperOutletEntropyMHD_hh
