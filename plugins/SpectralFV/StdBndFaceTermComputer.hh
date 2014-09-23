#ifndef COOLFluiD_Numerics_SpectralFV_StdBndFaceTermComputer_hh
#define COOLFluiD_Numerics_SpectralFV_StdBndFaceTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "SpectralFV/BaseBndFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a basic strategy that computes the boundary face terms
 * in a standard way using Gaussian quadrature
 *
 * @author Kris Van den Abeele
 */
class StdBndFaceTermComputer : public BaseBndFaceTermComputer {

public:  // methods

  /// Constructor
  StdBndFaceTermComputer(const std::string& name);

  /// Destructor
  ~StdBndFaceTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "StdBndFaceTermComputer";
  }

  std::string getType()
  {
    return "Std";
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
  void reconstructFaceAvgState(const std::vector< Framework::State* >& cellIntStates);

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

  /// quadrature wheights on SV faces
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > m_qWheightsSVFaces;

  /// polynomial values averaged over SV faces
  Common::SafePtr< std::vector< std::vector< CFreal > > > m_avgPolyValsSVFaces;

}; // class StdBndFaceTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFV_StdBndFaceTermComputer_hh
