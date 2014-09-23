#ifndef COOLFluiD_Numerics_SpectralFD_LESIPBndFaceTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_LESIPBndFaceTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "LES/LESVarSet.hh"

#include "SpectralFDNavierStokes/NSIPBndFaceTermComputer.hh"

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
class LESIPBndFaceTermComputer : public NSIPBndFaceTermComputer {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,LESIPBndFaceTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  LESIPBndFaceTermComputer(const std::string& name);

  /// Destructor
  ~LESIPBndFaceTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "LESIPBndFaceTermComputer";
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

}; // class LESIPBndFaceTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_LESIPBndFaceTermComputer_hh

