#ifndef COOLFluiD_Numerics_SpectralFD_LESIPFaceTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_LESIPFaceTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "LES/LESVarSet.hh"

#include "SpectralFDNavierStokes/NSIPFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the face terms for the Navier-Stokes/LES equations.
 *
 * @author Matteo Parsani
 * @author Kris Van den Abeele
 */
class LESIPFaceTermComputer : public NSIPFaceTermComputer {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,LESIPFaceTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  LESIPFaceTermComputer(const std::string& name);

  /// Destructor
  ~LESIPFaceTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "LESIPFaceTermComputer";
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
   * @pre cast m_faceDiffFluxComputer
   *      reconstructFluxPntsStates(), reconstructFluxPntsGradients(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  virtual void computeDiffFaceTerm(std::vector< RealVector >& resUpdates);


protected: // data

  /// Filter width volumes use in the LES calculation
  RealVector m_filterWidthVolumes;

}; // class LESIPFaceTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_LESIPFaceTermComputer_hh

