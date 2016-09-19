#ifndef COOLFluiD_Numerics_SpectralFV_BaseVolTermComputer_hh
#define COOLFluiD_Numerics_SpectralFV_BaseVolTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "SpectralFV/ReconstructStatesSpectralFV.hh"
#include "SpectralFV/SpectralFVElementData.hh"
#include "SpectralFV/SpectralFVMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a basic strategy that computes the volume terms in a cell
 *
 * @author Kris Van den Abeele
 */
class BaseVolTermComputer : public SpectralFVMethodStrategy {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFVMethodData,BaseVolTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  BaseVolTermComputer(const std::string& name);

  /// Destructor
  ~BaseVolTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BaseVolTermComputer";
  }
  
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// set the current cell
  void setCurrentCell(Framework::GeometricEntity* cell)
  {
    m_cell = cell;
  }

  /// Set up private data and data
  virtual void setup();

  /**
   * set volume term data for current element type
   */
  virtual void setVolumeTermData(CFuint iElemType) = 0;

  /**
   * compute cell data
   * @pre setCurrentCell
   */
  void computeCellData();

  /**
   * reconstruct solution in the required points
   * @pre setVolumeTermData()
   */
  void reconstructStates(const std::vector< Framework::State* >& cellStates);

  /**
   * reconstruct gradients in the required points
   * @pre setVolumeTermData()
   */
  void reconstructGradients(const std::vector< std::vector< RealVector >* >& cellGradients,
                                    const CFuint nbrCellGrads);

  /**
   * backup and reconstruct physical variable in the required points
   * @pre setVolumeTermData()
   */
  void backupAndReconstructPhysVar(const CFuint iVar, const std::vector< Framework::State* >& cellStates);

  /**
   * backup physical variable in the required points
   */
  void backupPhysVar(const CFuint iVar);

  /**
   * restore physical variable in the required points
   */
  void restorePhysVar(const CFuint iVar);

  /**
   * compute the convective volume term for this cell
   * @pre reconstructStates(), setVolumeTermData() and setCellFaceNormTransfM()
   */
  void computeCellConvVolumeTerm(RealVector& resUpdates);

  /**
   * volume term contribution to the gradients
   * @pre reconstructStates(), setVolumeTermData() and setCellFaceNormTransfM()
   */
  virtual void computeGradientVolumeTerm(std::vector< std::vector< RealVector > >& gradUpdates);

  /**
   * compute the diffusive volume term for this cell
   * @pre reconstructFaceStates(), setVolumeTermData(), setCellFaceNormTransfM() and
   *      setCellFaceTermData()
   */
  void computeCellDiffVolumeTerm(RealVector& resUpdates);

protected: // functions

  /**
   * compute the convective volume term from the solution in the flux points
   * @pre updateVarSet->computeStatesData
   */
  virtual void computeConvVolTermFromFlxPntSol(RealVector& resUpdates) = 0;

  /**
   * compute volume term contribution to the gradients from the solution in the flux points
   * @pre diffusiveVarSet->setGradientsVars()
   */
  virtual void computeGradVolTermFromFlxPntSol(std::vector< std::vector< RealVector > >& gradUpdates) = 0;

  /**
   * compute the diffusive volume term from the solution and the gradients in the flux points
   */
  virtual void computeDiffVolTermFromFlxPntSolAndGrad(RealVector& resUpdates) = 0;

protected: // data

  /// socket for face normal transformation matrices
  /// (this only applies to elements with a linear transformation to a reference element!!!)
  Framework::DataSocketSink<RealMatrix > socket_faceNormTransfMatrices;

  /// update variable set
  Common::SafePtr< Framework::ConvectiveVarSet > m_updateVarSet;

  /// Diffusive variable set
  Common::SafePtr< Framework::DiffusiveVarSet > m_diffusiveVarSet;

  /// Strategy that reconstructs the states in a given number of nodes
  Common::SafePtr< ReconstructStatesSpectralFV > m_statesReconstr;

  /// reconstruction coefficients for the flux points
  Common::SafePtr< std::vector< std::vector< CFreal > > > m_flxPntsRecCoefs;

  /// current cell
  Framework::GeometricEntity* m_cell;

  /// transformation matrix for cell normals
  RealMatrix m_cellFaceNormTransfM;

  /// variable for solution in flux points
  std::vector< Framework::State* > m_solInFlxPnts;

  /// variable for extra variables in flux points
  std::vector< RealVector* > m_extraVarsInFlxPnts;

  /// variable for solution in flux points stored in RealVector
  std::vector< RealVector* > m_solRVInFlxPnts;

  /// gradient variables in the flux points
  RealMatrix m_gradVarsInFlxPnts;

  /// gradients in the flux points
  std::vector< std::vector< RealVector* > > m_gradInFlxPnts;

  /// backup for physical variable in flux points
  std::vector< CFreal > m_backupPhysVar;

  /// number of flux points for current element type
  CFuint m_nbrFlxPnts;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

  /// number of extra variables in the physical model
  CFuint m_nbrExtraVars;

}; // class BaseVolTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFV_BaseVolTermComputer_hh

