#ifndef COOLFluiD_Numerics_SpectralFD_BR2BndFaceTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_BR2BndFaceTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "SpectralFD/CompactBndFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a basic strategy that computes the face terms at the boundaries
 *  using the Bassi-Rebay II scheme.
 *
 * @author Kris Van den Abeele
 */
class BR2BndFaceTermComputer : public CompactBndFaceTermComputer {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,BR2BndFaceTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  BR2BndFaceTermComputer(const std::string& name);

  /// Destructor
  ~BR2BndFaceTermComputer();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Configures the method, by allocating the it's dynamic members.
   *
   * @param args missing documentation
   */
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BR2BndFaceTermComputer";
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
   * reconstruct solution polynomial gradients in face flux points
   * @pre reconstructFluxPntsStates()
   */
  void reconstructFluxPntsSolPolyGrads(const std::vector< Framework::State* >& cellIntStates);

  /**
   * reconstruct given solution polynomial gradient in face flux points in the given cell
   * @pre backupAndReconstructPhysVar()
   */
  void backupAndReconstrFluxPntsSolPolyGrad(const CFuint iVar,
                                            const std::vector< Framework::State* >& cellStates);

  /**
   * backup physical variable gradient in the boundary face in the required points
   */
  void backupGradPhysVar(const CFuint iVar);

  /**
   * restore physical variable gradient in the boundary face in the required points
   */
  void restorePhysVarGrad(const CFuint iVar);

protected: // data

  /// solution point mapped coordinates
  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoords;

  /// backup for physical variable in flux points
  std::vector< std::vector< RealVector > > m_backupIntStateGrad;

  /// backup for the ghost state gradients
  std::vector< std::vector< RealVector > > m_backupGhostStateGrad;

  /// Jacobian determinants at solution points
  std::valarray<CFreal> m_solJacobDet;

  /// parameter of the BR2 approach determining the amount of damping
  CFreal m_alpha;

  /// lifting operator for the gradients in the solution points
  std::vector< std::vector< RealVector > > m_liftingOperators;

  /// pointers to lifting operator for the gradients in the solution points
  std::vector< std::vector< RealVector >* > m_liftingOperatorsPtrs;

  /// lifting operator terms for the gradients at the face flux points
  std::vector< std::vector< RealVector* > > m_liftOperatorTerms;

  /// lifting operator for a single gradient in the solution points
  std::vector< RealVector > m_liftingOperator;

  /// lifting operator for a single gradient in the solution points
  std::vector< std::vector< RealVector > > m_backupLiftingOperatorTerms;

}; // class BR2BndFaceTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_BR2BndFaceTermComputer_hh

