#ifndef COOLFluiD_FluxReconstructionMethod_ConvBndCorrectionsRHSFluxReconstruction_hh
#define COOLFluiD_FluxReconstructionMethod_ConvBndCorrectionsRHSFluxReconstruction_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"
#include "FluxReconstructionMethod/BaseCorrectionFunction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command that computes contribution of the boundary faces for the
   * Flux Reconstruction schemes for convective terms to the RHS
   *
   * @author Ray Vandenhoeck
   * @author Alexander Papen
   *
   */
class ConvBndCorrectionsRHSFluxReconstruction : public FluxReconstructionSolverCom {

public:
  typedef Framework::BaseMethodCommandProvider<
      FluxReconstructionSolverData,ConvBndCorrectionsRHSFluxReconstruction > PROVIDER;

public:

  /**
   * Constructor
   */
  ConvBndCorrectionsRHSFluxReconstruction(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~ConvBndCorrectionsRHSFluxReconstruction();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
 
  /**
   * Unset up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * set the BC state computer
   */
  void setBcStateComputer(Common::SafePtr< BCStateComputer > bcStateComputer)
  {
    m_bcStateComputer = bcStateComputer;
  }

protected: // functions

  /**
   * Execute on the current TRS
   */
  virtual void executeOnTrs();
  
  /// compute the states and ghost states in the flx pnts
  virtual void computeFlxPntStates();
  
  /// compute the interface flux
  virtual void computeInterfaceFlxCorrection();
  
  /// compute the total correction
  void computeCorrection(std::vector< RealVector >& corrections);

  /// add the residual updates to the RHS
  void updateRHS();

  /// compute the bnd face corrections to the gradients
  virtual void computeGradientBndFaceCorrections();
  
  /// add the updates to the wave speed
  void updateWaveSpeed();
  
  /// set the bnd face data necessary to compute FI-FD
  virtual void setBndFaceData(CFuint faceID);
  
  /**
   * compute the wave speed updates for this face
   * @pre reconstructFluxPntsStates(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  void computeWaveSpeedUpdates(CFreal& waveSpeedUpd);

protected: // data

  /// storage of the rhs
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// socket for updateCoeff
  /// denominators of the coefficients for the update
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// socket for gradients
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradients;
  
  /// socket for gradientsAV
  Framework::DataSocketSink< std::vector< RealVector > > socket_gradientsAV;
  
  /// socket for size of projection vector in face flux points
  Framework::DataSocketSink<  std::vector< CFreal > > socket_faceJacobVecSizeFaceFlxPnts;
  
  /// update variable set
  Common::SafePtr< Framework::ConvectiveVarSet > m_updateVarSet;

  /// builder of faces
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> > m_faceBuilder;

  /// the BCStateComputer for this BC
  Common::SafePtr< BCStateComputer > m_bcStateComputer;

  /// variable for current face
  Framework::GeometricEntity* m_face;

  /// variable for current cell
  Framework::GeometricEntity* m_intCell;

  /// variable for current face orientation
  CFuint m_orient;
  
  /// number of dimensions in the physical model
  CFuint m_dim;

  /// the states in the neighbouring cell
  std::vector< Framework::State* >* m_cellStates;

  /// update for the wave speed in the neighbouring cell
  CFreal m_waveSpeedUpd;

  /// number of equations in the physical model
  CFuint m_nbrEqs;
  
  /// number of solution pnts in the cell
  CFuint m_nbrSolPnts;
  
  /// number of flux pnts on a face
  CFuint m_nbrFaceFlxPnts;
  
  /// solution point mapped coordinates
  Common::SafePtr< std::vector< RealVector > > m_solPntsLocalCoords;
  
  // flx pnt mapped coordinates
  std::vector< RealVector > m_flxPntsLocalCoords;
  
  // All flux points of a cell
  Common::SafePtr<std::vector< RealVector > > m_allCellFlxPnts;
  
  /// flux point coordinates
  std::vector< RealVector > m_flxPntCoords;
  
  /// flx pnt - face connectivity
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_faceFlxPntConn;
  
  /// face connectivity per orient
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_faceConnPerOrient;
  
  /// coefficients for integration over a face
  Common::SafePtr< RealVector > m_faceIntegrationCoefs;
  
  /// face Jacobian vector sizes (abs)
  std::vector< CFreal > m_faceJacobVecAbsSizeFlxPnts;
  
  /// extrapolated states in the flux points of the cell
  std::vector< Framework::State* > m_cellStatesFlxPnt;
  
  /// ghost flux point solutions
  std::vector< Framework::State* > m_flxPntGhostSol;
  
  /// Riemann flux
  Common::SafePtr< RiemannFlux > m_riemannFluxComputer;
  
  /// Interface fluxes at the flux points of a face
  std::vector< RealVector> m_flxPntRiemannFlux;
  
  /// Correction function computer
  Common::SafePtr< BaseCorrectionFunction > m_corrFctComputer;
  
  /// Divergence of the correction function
  std::vector< std::vector< CFreal > > m_corrFctDiv;
  
  /// corrections to be added to RHS because of the boundary faces
  std::vector< RealVector> m_corrections;
  
  /// local cell face - mapped coordinate direction per orientation
  Common::SafePtr< std::vector< CFint > > m_faceMappedCoordDir;
  
  /// unit normal vector in flux points
  std::vector< RealVector > m_unitNormalFlxPnts;
  
  /// face Jacobian vector sizes
  std::vector< CFreal > m_faceJacobVecSizeFlxPnts;
  
  /// updates to the gradients
  std::vector< std::vector< RealVector > > m_gradUpdates;
  
  /// coefs to extrapolate the states to the flx pnts
  Common::SafePtr< std::vector< std::vector< CFreal > > > m_solPolyValsAtFlxPnts;
  
  /// local coordinates of the flux points on one face
  Common::SafePtr< std::vector< RealVector > > m_flxLocalCoords;

  /// local coordinates of the flux points on one face per face type
  Common::SafePtr<std::vector< std::vector< RealVector > > > m_faceFlxPntsLocalCoordsPerType;

  /// dependencies of flx pnts on sol pnts
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_flxSolDep;

  /// nbr of sol pnts on which a flx pnt is dependent
  CFuint m_nbrSolDep;
  
  /// vector for the face jacobian vectors
  std::vector< RealVector > m_faceJacobVecs;
  
  /// correction term projected to a normal
  RealVector m_projectedCorr;
  
  /// FR order
  CFuint m_order;
  
  /// Factor correcting Face normals direction (-1 factor needed for Tetra, due to the numbering convention the face normals are pointing inwards)
  CFreal m_mappedFaceNormalDir;

  private:

  /// Physical data temporary vector
  RealVector m_pData;

}; // end of class ConvBndCorrectionsRHSFluxReconstruction

//////////////////////////////////////////////////////////////////////////////

 } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ConvBndCorrectionsRHSFluxReconstruction_hh
