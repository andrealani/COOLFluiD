#ifndef COOLFluiD_Numerics_SpectralFD_IPBndFaceTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_IPBndFaceTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "SpectralFD/CompactBndFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a basic strategy that computes the face terms at the boundaries
 * if the approach is fully compact.
 *
 * @author Kris Van den Abeele
 */
class IPBndFaceTermComputer : public CompactBndFaceTermComputer {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,IPBndFaceTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  IPBndFaceTermComputer(const std::string& name);

  /// Destructor
  ~IPBndFaceTermComputer();

  /**
   * Configures the method, by allocating the it's dynamic members.
   *
   * @param args missing documentation
   */
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets the Class name
  static std::string getClassName()
  {
    return "IPBndFaceTermComputer";
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
   * backup physical variable gradient in the boundary face in the required points
   */
  void backupGradPhysVar(const CFuint iVar);

  /**
   * restore physical variable gradient in the boundary face in the required points
   */
  void restorePhysVarGrad(const CFuint iVar);

protected: // data

  /// backup for physical variable in flux points
  std::vector< RealVector > m_backupPhysVarGrad;

  /// backup for the ghost state gradients
  std::vector< std::vector< RealVector > > m_backupGhostStateGrad;

}; // class IPBndFaceTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_IPBndFaceTermComputer_hh

