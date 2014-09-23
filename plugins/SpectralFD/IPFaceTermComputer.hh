#ifndef COOLFluiD_Numerics_SpectralFD_IPFaceTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_IPFaceTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/MatrixInverter.hh"

#include "SpectralFD/CompactFaceTermComputer.hh"
#include "SpectralFD/FaceDiffusiveFluxIPApproach.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a basic strategy that computes the face terms
  * with the interior penalty approach.
 *
 * @author Kris Van den Abeele
 */
class IPFaceTermComputer : public CompactFaceTermComputer {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,IPFaceTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  IPFaceTermComputer(const std::string& name);

  /// Destructor
  ~IPFaceTermComputer();

  /**
   * Configures the method, by allocating the it's dynamic members.
   *
   * @param args missing documentation
   */
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets the Class name
  static std::string getClassName()
  {
    return "IPFaceTermComputer";
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
   * compute neighbour cell data
   * @pre setCurrentFace
   */
  virtual void computeNeighbourCellData();

  /**
   * backup physical variable gradient in one cell in the required points
   */
  void backupPhysVarGrad(const CFuint side, const CFuint iVar);

  /**
   * restore physical variable gradient in one cell in the required points
   */
  void restorePhysVarGrad(const CFuint side, const CFuint iVar);

protected: // data

  /// backup for state gradients
  std::vector< RealVector > m_backupPhysVarGrad;

}; // class IPFaceTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_IPFaceTermComputer_hh

