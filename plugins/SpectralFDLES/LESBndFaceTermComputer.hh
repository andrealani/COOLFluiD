#ifndef COOLFluiD_Numerics_SpectralFD_LESBndFaceTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_LESBndFaceTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "LES/LESVarSet.hh"

#include "SpectralFDNavierStokes/NavierStokesBndFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the face terms at the boundaries for the Navier-Stokes/LES equations.
 *
 * @author Kris Van den Abeele
 * @author Matteo Parsani
 */
class LESBndFaceTermComputer : public NavierStokesBndFaceTermComputer {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,LESBndFaceTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  LESBndFaceTermComputer(const std::string& name);

  /// Destructor
  ~LESBndFaceTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "LESBndFaceTermComputer";
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

  /// LES variable set
  Common::SafePtr< LES::LESVarSet > m_lesVarSet;

  /// Filter width volumes to use in the LES calculation
  RealVector m_filterWidthVolumes;

}; // class LESBndFaceTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_LESBndFaceTermComputer_hh

