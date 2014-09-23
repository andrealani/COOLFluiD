#ifndef COOLFluiD_Numerics_SpectralFD_LESFaceTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_LESFaceTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "LES/LESVarSet.hh"

#include "SpectralFDNavierStokes/NavierStokesFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a strategy that computes the face terms for the Navier-Stokes/LES equations.
 *
 * @author Kris Van den Abeele
 * @author Matteo Parsani
 */
class LESFaceTermComputer : public NavierStokesFaceTermComputer {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,LESFaceTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  LESFaceTermComputer(const std::string& name);

  /// Destructor
  ~LESFaceTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "LESFaceTermComputer";
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

  /// LES variable set
  Common::SafePtr< LES::LESVarSet > m_lesVarSet;

  /// Filter width volumes use in the LES calculation
  RealVector m_filterWidthVolumes;

}; // class LESFaceTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_LESFaceTermComputer_hh

