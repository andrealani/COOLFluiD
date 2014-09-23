#ifndef COOLFluiD_Numerics_SpectralFD_LESVolTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_LESVolTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "LES/LESVarSet.hh"

#include "SpectralFDNavierStokes/NavierStokesVolTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the volume terms in a cell for the Navier-Stokes/LES equations
 *
 * @author Kris Van den Abeele
 * @author Matteo Parsani
 */
class LESVolTermComputer : public NavierStokesVolTermComputer {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,LESVolTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  LESVolTermComputer(const std::string& name);

  /// Destructor
  ~LESVolTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "LESVolTermComputer";
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

}; // class LESVolTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_LESVolTermComputer_hh

