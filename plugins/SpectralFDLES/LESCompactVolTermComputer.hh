#ifndef COOLFluiD_Numerics_SpectralFD_LESCompactVolTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_LESCompactVolTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "LES/LESVarSet.hh"

#include "SpectralFDNavierStokes/NSCompactVolTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the volume terms in a cell for the Navier-Stokes/LES equations
 *
 * @author Matteo Parsani
 * @author Kris Van den Abeele
 */
class LESCompactVolTermComputer : public NSCompactVolTermComputer {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,LESCompactVolTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  LESCompactVolTermComputer(const std::string& name);

  /// Destructor
  ~LESCompactVolTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "LESCompactVolTermComputer";
  }

  /// Setup private data
  virtual void setup();

  /// Unsetup private data
  virtual void unsetup();

  /**
   * set volume term data for current element type
   */
  virtual void setVolumeTermData(CFuint iElemType);

  /**
   * compute cell data
   * @pre setCurrentCell
   */
  virtual void computeCellData();

protected: // functions
 /**
   * compute the diffusive volume term for this cell from the solution and the gradients in the flux points
   * @pre reconstructFaceStates(), setVolumeTermData(), setCellFaceNormTransfM() and
   *      setCellFaceTermData()
   */
void computeDiffVolTermFromFlxPntSolAndGrad(RealVector& resUpdates);

protected: // data

  /// LES variable set
  Common::SafePtr< LES::LESVarSet > m_lesVarSet;

  /// Number of internal flux points
  CFreal m_nbrOfIntFlxPnts;

  /// Filter width volumes to use in the LES calculation
  RealVector m_filterWidthVolumes;

}; // class LESCompactVolTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_LESCompactVolTermComputer_hh

