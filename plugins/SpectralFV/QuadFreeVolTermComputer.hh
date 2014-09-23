#ifndef COOLFluiD_Numerics_SpectralFV_QuadFreeVolTermComputer_hh
#define COOLFluiD_Numerics_SpectralFV_QuadFreeVolTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "SpectralFV/BaseVolTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the volume terms in a cell
 * using a quadrature-free approach.
 *
 * @author Kris Van den Abeele
 */
class QuadFreeVolTermComputer : public BaseVolTermComputer {

public:  // methods

  /// Constructor
  QuadFreeVolTermComputer(const std::string& name);

  /// Destructor
  ~QuadFreeVolTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "QuadFreeVolTermComputer";
  }

  /// Set up private data and data
  void setup();

  /**
   * Unsetup private data
   */
  void unsetup();

  /**
   * set volume term data for current element type
   */
  void setVolumeTermData(CFuint iElemType);

protected: // functions

  /**
   * compute the convective volume term from the solution in the flux points
   * @pre updateVarSet->computeStatesData
   */
  void computeConvVolTermFromFlxPntSol(RealVector& resUpdates);

  /**
   * compute volume term contribution to the gradients from the solution in the flux points
   * @pre diffusiveVarSet->setGradientsVars()
   */
  void computeGradVolTermFromFlxPntSol(std::vector< std::vector< RealVector > >& gradUpdates);

  /**
   * compute the diffusive volume term from the solution and the gradients in the flux points
   */
  void computeDiffVolTermFromFlxPntSolAndGrad(RealVector& resUpdates);

protected: // data

  /// tensor for the evaluation of the volume terms
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > m_volTermTensor;

  /// physical flux in a flux point
  RealMatrix m_physFlux;

  /// variable for flux through a face
  RealVector m_faceFlux;

  /// number of equations in the physical model
  CFuint m_dim;

}; // class QuadFreeVolTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFV_QuadFreeVolTermComputer_hh

