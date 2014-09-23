#ifndef COOLFluiD_Numerics_SpectralFV_QuadFreeFaceTermComputer_hh
#define COOLFluiD_Numerics_SpectralFV_QuadFreeFaceTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "SpectralFV/BaseFaceTermComputer.hh"
#include "SpectralFV/SpectralFVMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a basic strategy that computes the face terms
 * in using a quadrature free method
 *
 * @author Kris Van den Abeele
 */
class QuadFreeFaceTermComputer : public BaseFaceTermComputer {

public:  // methods

  /// Constructor
  QuadFreeFaceTermComputer(const std::string& name);

  /// Destructor
  ~QuadFreeFaceTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "QuadFreeFaceTermComputer";
  }

  /// Set up private data and data
  void setup();

  /// Unset up private data and data
  void unsetup();

  /**
   * set face term data
   */
  void setFaceTermData();

  /**
   * reconstruct face averaged solutions
   * @pre setFaceTermData() and reconstructFluxPntsStates()
   */
  void reconstructFaceAvgState(const std::vector< Framework::State* >& cellLStates,
                                       const std::vector< Framework::State* >& cellRStates);

protected: // functions

  /**
   * compute the flux integrals from the fluxes in the flux points
   * @pre fluxes in flux points are computed
   */
  void computeFaceFluxIntegralFromFlxPntFluxes(RealVector& resUpdates);

  /**
   * compute face term contribution to the gradients from the solution in the flux points
   * @pre diffusiveVarSet->setGradientsVars()
   */
  void computeGradFaceTermFromFlxPntSol(std::vector< std::vector< RealVector > >& gradUpdates);

protected: // data

  /// coefficients for the computation of the flux through an external CV face
  Common::SafePtr< std::vector< std::vector< CFreal > > > m_cvExtFaceFluxCoef;

  /// coefficients for the computation of the average solution over an external CV face
  Common::SafePtr< std::vector< CFreal > > m_avgSolInSVFaceCoef;

}; // class QuadFreeFaceTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFV_QuadFreeFaceTermComputer_hh
