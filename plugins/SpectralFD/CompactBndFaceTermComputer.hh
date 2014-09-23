#ifndef COOLFluiD_Numerics_SpectralFD_CompactBndFaceTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_CompactBndFaceTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "MathTools/MatrixInverter.hh"

#include "SpectralFD/BaseBndFaceTermComputer.hh"

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
class CompactBndFaceTermComputer : public BaseBndFaceTermComputer {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,CompactBndFaceTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  CompactBndFaceTermComputer(const std::string& name);

  /// Destructor
  ~CompactBndFaceTermComputer();

  /**
   * Configures the method, by allocating the it's dynamic members.
   *
   * @param args missing documentation
   */
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets the Class name
  static std::string getClassName()
  {
    return "CompactBndFaceTermComputer";
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
  virtual void reconstructFluxPntsSolPolyGrads(const std::vector< Framework::State* >& cellIntStates);

  /**
   * reconstruct given solution polynomial gradient in face flux points in the given cell
   */
  virtual void backupAndReconstrFluxPntsSolPolyGrad(const CFuint iVar,
                                                    const std::vector< Framework::State* >& cellStates);

  /**
   * backup physical variable gradient in the boundary face in the required points
   */
  virtual void backupGradPhysVar(const CFuint iVar) = 0;

  /**
   * restore physical variable gradient in the boundary face in the required points
   */
  virtual void restorePhysVarGrad(const CFuint iVar) = 0;

  /**
   * set cell mapped coordinates of point set in which the states should be reconstructed
   * @pre setFaceOrientation()
   */
  virtual void setPointSet(const std::vector< RealVector >& faceMappedCoords);

  /**
   * compute face data for given set of points
   * @pre setCurrentFace and setOutputData()
   */
  virtual void computeFacePntSetData();

  /**
   * reconstruct gradients in point set
   * @pre setOutputData()
   */
  virtual std::vector< std::vector< RealVector > >& reconstructGivenPntsGrads(const std::vector< Framework::State* >& cellIntStates);

protected: // data

  /// flux point solution polynomial derivative coefficients
  Common::SafePtr< std::vector< std::vector< std::vector< CFreal > > > > m_coefSolPolyDerivInFlxPnts;

  /// flux point solution polynomial derivative indexes
  Common::SafePtr< std::vector< std::vector< std::vector< CFuint > > > > m_solPntIdxsSolPolyDerivInFlxPnts;

  /// cell inverse Jacobian matrices
  std::vector< RealMatrix > m_cellInvJacobMatr;

  /// helper matrix
  RealMatrix m_auxMatrix;

  /// inverter of matrices with size dim x dim
  MathTools::MatrixInverter* m_matrInverter;

  /// solution polynomial derivation coefficients of a given set of points
  std::vector< std::vector< std::vector< CFreal > > > m_solPolyDerivCoefsPntSet;

  /// cell inverse Jacobian matrices in given point set
  std::vector< RealMatrix > m_cellInvJacobMatrPntSet;

}; // class CompactBndFaceTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_CompactBndFaceTermComputer_hh

