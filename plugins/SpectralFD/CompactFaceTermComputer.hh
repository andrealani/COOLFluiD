#ifndef COOLFluiD_Numerics_SpectralFD_CompactFaceTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_CompactFaceTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/MatrixInverter.hh"

#include "SpectralFD/BaseFaceTermComputer.hh"
#include "SpectralFD/FaceDiffusiveFluxCompactApproach.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a basic strategy that computes the face terms if the scheme is fully compact
 *
 * @author Kris Van den Abeele
 */
class CompactFaceTermComputer : public BaseFaceTermComputer {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,CompactFaceTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  CompactFaceTermComputer(const std::string& name);

  /// Destructor
  ~CompactFaceTermComputer();

  /**
   * Configures the method, by allocating the it's dynamic members.
   *
   * @param args missing documentation
   */
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets the Class name
  static std::string getClassName()
  {
    return "CompactFaceTermComputer";
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
   */
  virtual void reconstructFluxPntsSolPolyGrads(const std::vector< std::vector< Framework::State* >* >& cellStates);

  /**
   * reconstruct given solution polynomial gradient in face flux points in the given cell
   */
  virtual void backupAndReconstrFluxPntsSolPolyGrad(const CFuint side,
                                                    const CFuint iVar,
                                                    const std::vector< Framework::State* >& cellStates);

  /**
   * backup physical variable gradient in one cell in the required points
   */
  virtual void backupPhysVarGrad(const CFuint side, const CFuint iVar) = 0;

  /**
   * restore physical variable gradient in one cell in the required points
   */
  virtual void restorePhysVarGrad(const CFuint side, const CFuint iVar) = 0;

  /**
   * face term contribution to a physical variable gradient
   * @pre reconstructFluxPntsStates(), setFaceTermData() and set the geometrical data of the face
   */
  void computePhysVarGradFaceTerm(const CFuint iVar, std::vector< std::vector< RealVector > >& physVarGradUpdates);

protected: // functions

  /**
   * compute face term contribution to a physical variable gradient from the solution in the flux points
   */
  void computePhysVarGradFaceTermFromFlxPntSol(const CFuint iVar, std::vector< std::vector< RealVector > >& physVarGradUpdates);

protected: // data

  /// Face diffusive flux
  Common::SafePtr< FaceDiffusiveFluxCompactApproach > m_faceCompactDiffFluxComputer;

  /// flux point solution polynomial derivative coefficients
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > m_coefSolPolyDerivInFlxPnts;

  /// flux point solution polynomial derivative indexes
  Common::SafePtr< std::vector< std::vector< std::vector< CFuint > > > > m_solPntIdxsSolPolyDerivInFlxPnts;

  /// cell inverse Jacobian matrices
  std::vector< std::vector< RealMatrix > > m_cellInvJacobMatr;

  /// helper matrix
  RealMatrix m_auxMatrix;

  /// inverter of matrices with size dim x dim
  MathTools::MatrixInverter* m_matrInverter;

}; // class CompactFaceTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_CompactFaceTermComputer_hh

