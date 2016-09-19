#ifndef COOLFluiD_Numerics_SpectralFD_BaseVolTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_BaseVolTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "SpectralFD/ReconstructStatesSpectralFD.hh"
#include "SpectralFD/SpectralFDMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a basic strategy that computes the volume terms in a cell
 *
 * @author Kris Van den Abeele
 */
class BaseVolTermComputer : public SpectralFDMethodStrategy {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,BaseVolTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  BaseVolTermComputer(const std::string& name);

  /// Destructor
  ~BaseVolTermComputer();

  /**
   * Configures the method, by allocating the it's dynamic members.
   *
   * @param args missing documentation
   */
  virtual void configure ( Config::ConfigArgs& args );

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BaseVolTermComputer";
  }
  
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
  /// set the current cell
  void setCurrentCell(Framework::GeometricEntity* cell)
  {
    m_cell = cell;
  }

  /// @return m_cellExtraVars
  Common::SafePtr< std::vector< RealVector* > > getCellExtraVars()
  {
    return &m_cellExtraVars;
  }

  /// Setup private data
  virtual void setup();

  /// Unsetup private data
  virtual void unsetup();

  /**
   * set volume term data for current element type
   */
  virtual void setVolumeTermData(CFuint iElemType);

  /**
   * compute cell data
   * @pre setCurrentCell
   */
  virtual void computeCellData();

  /**
   * reconstruct solution in the required points
   * @pre setVolumeTermData()
   */
  void reconstructStates(const std::vector< Framework::State* >& cellStates, bool onlyExtraVars = false);

  /**
   * reconstruct gradients in the required points
   * @pre setVolumeTermData()
   */
  void reconstructGradients(const std::vector< std::vector< RealVector >* >& cellGradients);

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
  virtual void computeCellDiffVolumeTerm(RealVector& resUpdates);

  /**
   * volume term contribution to the gradients of the extra variables
   * @pre reconstructStates(), setVolumeTermData() and setCellFaceNormTransfM()
   */
  virtual void computeGradientExtraVarsVolumeTerm(std::vector< std::vector< RealVector > >& gradUpdates);

  /**
   * Returns the DataSocket's that this strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // functions

  /**
   * compute the convective volume term from the solution in the flux points
   * @pre updateVarSet->computeStatesData
   */
  void computeConvVolTermFromFlxPntSol(RealVector& resUpdates);

  /**
   * compute volume term contribution to the gradients from the solution in the flux points
   * @pre diffusiveVarSet->setGradientsVars()
   */
  void computeGradVolTermFromFlxPntSol(std::vector< std::vector< RealVector > >& gradUpdates);

  /**
   * compute the diffusive volume term from the solution and the gradients in the flux points
   */
  virtual void computeDiffVolTermFromFlxPntSolAndGrad(RealVector& resUpdates);

protected: // data

  /// the socket containing the extra variables
  Framework::DataSocketSink< RealVector > socket_extraVars;

  /// update variable set
  Common::SafePtr< Framework::ConvectiveVarSet > m_updateVarSet;

  /// Diffusive variable set
  Common::SafePtr< Framework::DiffusiveVarSet > m_diffusiveVarSet;

  /// Strategy that reconstructs the states in a given number of nodes
  Common::SafePtr< ReconstructStatesSpectralFD > m_statesReconstr;

  /// reconstruction coefficients for the flux points
  Common::SafePtr< RealMatrix > m_flxPntsRecCoefs;

  /// derivation coefficients for the solution points
  Common::SafePtr< RealMatrix > m_solPntsDerivCoefs;

  /// indexes of internal flux points
  Common::SafePtr< std::vector< CFuint > > m_intFlxPntIdxs;

  /// derivation direction in the flux points
  Common::SafePtr< std::vector< CFuint > > m_flxPntDerivDir;

  /// derivation direction in the internal flux points
  std::vector< CFuint > m_intFlxPntDerivDir;

  /// flux points mapped coordinates
  Common::SafePtr< std::vector< RealVector > > m_flxPntsLocalCoords;

  /// mapped coordinates of the internal flux points
  std::vector< RealVector > m_intFlxPntMappedCoord;

  /// flux point index (in the matrix flxPntRecCoefs) for reconstruction
  Common::SafePtr< std::vector< CFuint > > m_flxPntMatrixIdxForReconstruction;

  /// solution point index (in the cell) for reconstruction
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_solPntIdxsForReconstruction;

  /// flux point index (in the matrix m_solPntsDerivCoefs) for derivation
  Common::SafePtr< std::vector< CFuint > > m_flxPntMatrixIdxForDerivation;

  /// solution point index (in the cell) for derivation
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_solPntIdxsForDerivation;

  /// current cell
  Framework::GeometricEntity* m_cell;

  /// flux projection vectors in flux points
  std::vector< RealVector > m_cellFluxProjVects;

  /// cell extra variables
  std::vector< RealVector* > m_cellExtraVars;

  /// variable for solution in flux points
  std::vector< Framework::State* > m_solInFlxPnts;

  /// variable for solution in flux points stored in RealVector
  std::vector< RealVector* > m_solRVInFlxPnts;

  /// variable for extra variables in flux points
  std::vector< RealVector* > m_extraVarsInFlxPnts;

  /// gradient variables in the flux points
  RealMatrix m_gradVarsInFlxPnts;

  /// gradients in the flux points
  std::vector< std::vector< RealVector* > > m_gradInFlxPnts;

  /// gradient variable gradients in flux points
  std::vector< std::vector< RealVector* > > m_gradVarGradsInFlxPnts;

  /// backup for physical variable in flux points
  std::vector< CFreal > m_backupPhysVar;

  /// number of flux points for current element type
  CFuint m_nbrFlxPnts;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

  /// dimensionality
  CFuint m_dim;

  /// helper term for gradient computation
  RealVector m_gradTerm;

  /// number of extra variables in the physical model
  CFuint m_nbrExtraVars;
  
private:

  /// Physical data temporary vector
  RealVector m_pData;

}; // class BaseVolTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_BaseVolTermComputer_hh

