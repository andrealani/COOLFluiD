#ifndef COOLFluiD_Numerics_SpectralFD_LESBR2BndFaceTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_LESBR2BndFaceTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "LES/LESVarSet.hh"

#include "SpectralFDNavierStokes/NSBR2BndFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the face terms at the boundaries for the Navier-Stokes/LES equations.
 *
 * @author Matteo Parsani
 * @author Kris Van den Abeele
 */
class LESBR2BndFaceTermComputer : public NSBR2BndFaceTermComputer {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,LESBR2BndFaceTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  LESBR2BndFaceTermComputer(const std::string& name);

  /// Destructor
  ~LESBR2BndFaceTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "LESBR2BndFaceTermComputer";
  }

  /**
   * Set up private data and data
   */
  virtual void setup();

  /**
   * Unset up private data and data
   */
  virtual void unsetup();

  /**
   * set face term data
   */
  virtual void setFaceTermData();

  /**
   * compute neighbour cell data
   * @pre setCurrentFace
   */
  virtual void computeNeighbourCellData();

   /**
   * compute the diffusive face term for this face
   * @pre reconstructFluxPntsStates(), reconstructFluxPntsGradients(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  virtual void computeDiffFaceTerm(RealVector& resUpdates);

protected: // data

  /// Filter width volumes to use in the LES calculation
  RealVector m_filterWidthVolumes;

}; // class LESBR2BndFaceTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_LESBR2BndFaceTermComputer_hh

