#ifndef COOLFluiD_Numerics_SpectralFD_BR2FaceTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_BR2FaceTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/CompactFaceTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a basic strategy that computes the face terms using the Bassi-Rebay II scheme.
 *
 * @author Kris Van den Abeele
 */
class BR2FaceTermComputer : public CompactFaceTermComputer {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,BR2FaceTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  BR2FaceTermComputer(const std::string& name);

  /// Destructor
  ~BR2FaceTermComputer();

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
    return "BR2FaceTermComputer";
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
  void reconstructFluxPntsSolPolyGrads(const std::vector< std::vector< Framework::State* >* >& cellStates);

  /**
   * reconstruct given solution polynomial gradient in face flux points in the given cell
   * @pre backupAndReconstructPhysVar()
   */
  void backupAndReconstrFluxPntsSolPolyGrad(const CFuint side,
                                            const CFuint iVar,
                                            const std::vector< Framework::State* >& cellStates);

  /**
   * backup physical variable gradient in one cell in the required points
   */
  void backupPhysVarGrad(const CFuint side, const CFuint iVar);

  /**
   * restore physical variable gradient in one cell in the required points
   */
  void restorePhysVarGrad(const CFuint side, const CFuint iVar);

protected: // data

  /// solution point mapped coordinates
  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoords;

  /// Jacobian determinants at solution points
  std::vector< std::valarray<CFreal> > m_solJacobDet;

  /// backup for state gradients
  std::vector< std::vector< RealVector > > m_backupPhysVarGrad;

  /// parameter of the BR2 approach determining the amount of damping
  CFreal m_alpha;

  /// lifting operator for the gradients in the solution points
  std::vector< std::vector< std::vector< RealVector > > > m_liftingOperators;

  /// pointers to lifting operator for the gradients in the solution points
  std::vector< std::vector< std::vector< RealVector >* > > m_liftingOperatorsPtrs;

  /// lifting operator terms for the gradients at the face flux points
  std::vector< std::vector< std::vector< RealVector* > > > m_liftOperatorTerms;

  /// lifting operator for a single gradient in the solution points
  std::vector< std::vector< RealVector > > m_liftingOperator;

  /// lifting operator for a single gradient in the solution points
  std::vector< std::vector< RealVector > > m_backupLiftingOperatorTerms;

}; // class BR2FaceTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_BR2FaceTermComputer_hh

