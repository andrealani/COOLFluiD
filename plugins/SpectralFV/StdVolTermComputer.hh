#ifndef COOLFluiD_Numerics_SpectralFV_StdVolTermComputer_hh
#define COOLFluiD_Numerics_SpectralFV_StdVolTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"


#include "SpectralFV/BaseVolTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the volume terms in a cell
 * in a standard way using Gaussian quadrature.
 *
 * @author Kris Van den Abeele
 */
class StdVolTermComputer : public BaseVolTermComputer {

public:  // methods

  /// Constructor
  StdVolTermComputer(const std::string& name);

  /// Destructor
  ~StdVolTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "StdVolTermComputer";
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

  /// local internal face-CV connectivity
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_localIntFaceCVConn;

  /// local internal face normals in a SV
  Common::SafePtr< std::vector< std::vector< RealVector > > > m_intFaceQuadPntNorm;

  /// local internal face quadrature wheights
  Common::SafePtr< std::vector< std::vector< CFreal > > > m_intFaceQuadWheights;

  /// variable for normal vectors
  RealVector m_normal;

  /// variable for flux through a face
  RealVector m_faceFlux;

  /// helper term for gradient computation
  RealVector m_gradTerm;

}; // class StdVolTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFV_StdVolTermComputer_hh

