// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_ConvDiffLLAVJacobFluxReconstruction_hh
#define COOLFluiD_FluxReconstructionMethod_ConvDiffLLAVJacobFluxReconstruction_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"
#include "FluxReconstructionMethod/BaseCorrectionFunction.hh"
#include "FluxReconstructionMethod/DiffRHSJacobFluxReconstruction.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// This is a standard command to assemble the convective, diffusive 
/// and artificial viscosity part of the system using a FluxReconstruction solver for an implicit scheme
/// @author Ray Vandenhoeck
class ConvDiffLLAVJacobFluxReconstruction : public DiffRHSJacobFluxReconstruction {

public: // functions

  /// Constructor
  explicit ConvDiffLLAVJacobFluxReconstruction(const std::string& name);

  /// Destructor
  virtual ~ConvDiffLLAVJacobFluxReconstruction() {}

  /// Execute processing actions
  virtual void execute();
  
  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unsetup private data
   */
  virtual void unsetup();
  
  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();
  
  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
    needsSockets();
    
protected: //functions
    
  /// compute the interface flux
  virtual void computeInterfaceFlxCorrection();
  
  /**
   * Set the data for the current face necessary to calculate FI
   */
  virtual void setFaceData(CFuint faceID);
  
  /**
   * Compute the left and right states and gradients in the flx pnts
   */
  virtual void computeFlxPntStatesAndGrads();
  
  /**
   * compute the wave speed updates for this face
   * @pre reconstructFluxPntsStates(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  virtual void computeWaveSpeedUpdates(std::vector< CFreal >& waveSpeedUpd);
  
  /**
   * compute the unperturbed cell diffusive residuals
   * @pre m_faceTermComputers->computeDiffFaceTermAndUpdateCoefContributions
   * @pre setCellsData()
   */
  virtual void computeUnpertCellDiffResiduals(const CFuint side);

  /**
   * compute the contribution of the diffusive face term to both Jacobians
   */
  void computeBothJacobsDiffFaceTerm();

  /**
   * compute the contribution of the diffusive face term to one Jacobians
   */
  void computeOneJacobDiffFaceTerm(const CFuint side);
  
  /**
   * Compute the projected states on order P-1
   */
  void computeProjStates(std::vector< RealVector >& projStates);
  
  /**
   * Compute the projected states on order P-1 on a side of the current face
   */
  void computeProjStates(std::vector< RealVector >& projStates, const CFuint side);
  
  /**
   * Compute the artificial viscosity
   */
  virtual void computeEpsilon();
  
  /**
   * Compute the reference artificial viscosity
   */
  virtual void computeEpsilon0();
  
  /**
   * Compute the reference artificial viscosity
   */
  virtual void computeEpsilon0(const CFuint side);
  
  /**
   * Compute the smoothness indicator
   */
  virtual void computeSmoothness();
  
  /**
   * Compute the smoothness indicator
   */
  virtual void computeSmoothness(const CFuint side);
  
  /**
   * Store the computed artificial viscosity
   */
  virtual void storeEpsilon();
  
  /**
   * compute the artificial diffusive flux
   */
  virtual void computeFlux(const RealVector& values, const std::vector< RealVector* >& gradients, const RealVector& normal, const CFreal& radius, RealVector& flux);
  
  /**
   * command to compute the peclet number based on user input
   */
  virtual CFreal computePeclet(); 
  
  /// compute the volume term contribution to the gradients
  virtual void computeGradients();
  
  /// compute the face correction to the corrected gradients
  virtual void computeGradientFaceCorrections();
  
  /// set face data for gradient computation
  void setFaceDataForGradients(const CFuint faceID);
  
  /// initialize the data needed for the jacobian
  void initJacobianComputation();
  
  /// compute the cell flux jacobian numerically
  void computeCellFluxJacobianNum(const CFreal resFactor);
  
  /// compute the Riemann flux jacobian numerically
  virtual void computeRiemannFluxJacobianNum(const CFreal resFactor);
  
  /// compute the flux to gradient jacobian numerically
  void computeFluxToGradJacobianNum(const CFreal resFactor);
  
  /// compute the gradient to state jacobian analytically
  virtual void computeGradToStateJacobianAna();
  
  /// compute the Riemann flux to gradient jacobian numerically
  virtual void computeRiemannFluxToGradJacobianNum(const CFreal resFactor);
  
  /// compute the gradient variables to state jacobians numerically
  virtual void computeGradVarsToStateJacobianNum();
  
  /// compute the AV to state Jacobian
  virtual void computeEpsToStateJacobianAna() = 0;
  
  /// compute the cell LLAV Jacobian (to state) analytically
  virtual void computeLLAVCellFluxJacobianAna(const CFreal resFactor);
  
  /// compute the Riemann LLAV Jacobian (to state) analytically
  virtual void computeLLAVRiemannFluxJacobianAna(const CFreal resFactor);
  
protected: //data
    
  /// update variable set
  Common::SafePtr< Framework::ConvectiveVarSet > m_updateVarSet;
    
  /// order of the FR method
  CFuint m_order;
  
  /// transformation matrices to order P-1
  RealMatrix m_transformationMatrix;
  
  /// states projected on P-1
  std::vector< RealVector > m_statesPMinOne;
  
  /// artificial Viscosity
  CFreal m_epsilon;
  
  /// artificial Viscosity in the sol pnts
  std::vector< std::vector< CFreal > > m_solEpsilons;
  
  /// artificial Viscosity
  std::vector< std::vector< CFreal > > m_epsilonLR;
  
  /// reference artificial Viscosity
  CFreal m_epsilon0;
  
  /// reference smoothness
  CFreal m_s0;
  
  /// smoothness
  CFreal m_s;
  
  /// controlling parameter kappa
  CFreal m_kappa;
  
  /// show rate of LLAV info
  CFuint m_showrate;
  
  /// peclet number
  CFreal m_peclet;
  
  /// damping coefficient of recalculation of epsilon
  CFreal m_dampingCoeff;
  
  /// average artificial viscosities in the nodes
  RealVector m_nodeEpsilons;
  
  /// number of neighbors for each node
  RealVector m_nbNodeNeighbors;
  
  /// average artificial viscosities in the elements
  RealVector m_cellEpsilons;
  
  /// vector containing pointers to the nodes in a cell
  std::vector< std::vector< Framework::Node*  >* > m_cellNodes;
  
  /// number of corner nodes for current element type
  CFuint m_nbrCornerNodes;
  
  /// vector containing pointers to the nodes in a face
  std::vector< Framework::Node*  >* m_faceNodes;
  
  /// flag telling whether to compute the number of node neighbors
  bool m_flagComputeNbNghb;
  
  /// polynomial coefficients for reconstruction of the artificial viscosity at the flx pnts
  std::vector< std::vector< CFreal > > m_nodePolyValsAtFlxPnts;
  
  /// polynomial coefficients for reconstruction of the artificial viscosity at the sol pnts
  std::vector< std::vector< CFreal > > m_nodePolyValsAtSolPnts;
  
  /// cell node connectivity table
  Common::SafePtr< Framework::MeshData::ConnTable > m_cellNodesConn;
  
  /// current element index
  CFuint m_elemIdx;
  
  /// variable for faces
  const std::vector< Framework::GeometricEntity* >* m_facesCell;
  
  /// bool telling whether the jacobian is being computed
  bool m_jacob;
  
  /// residual after which the limiter is frozen
  CFreal m_freezeLimiterRes;
  
  /// iteration after which the limiter is frozen
  CFuint m_freezeLimiterIter;
  
  /// boolean telling whether to use max artificial viscosity wrt previous iteration
  bool m_useMax;
  
  /// total artificial viscosity added
  CFreal m_totalEps;
  
  /// total artificial viscosity added over all CPUs
  CFreal m_totalEpsGlobal;
  
  /// maximum smoothness in domain
  CFreal m_Smax;
  
  /// maximum smoothness in domain over all CPUs
  CFreal m_SmaxGlobal;
  
  /// boolean telling whether to add the contribution of the artificial flux to the update coefficients
  bool m_addUpdCoeff;
  
  /// value at which point positivity preservation is added
  CFreal m_minValue;
  
  /// index of the monitored variable for LLAV
  CFuint m_monitoredVar;
  
  /// subcell resolution
  CFreal m_subcellRes;
  
  /// index of the monitored physical variable for LLAV
  CFuint m_monitoredPhysVar;
  
  /// extra vector to store the unit normal vectors in the flx pnts
  std::vector< RealVector > m_unitNormalFlxPnts2;
  
  /// storage for the artificial viscosity
  Framework::DataSocketSource<CFreal> socket_artVisc;
  
  /// storage for the monitored phys var
  Framework::DataSocketSource<CFreal> socket_monPhysVar;
  
  /// storage for the smoothness
  Framework::DataSocketSource<CFreal> socket_smoothness;
  
  /// storage for the normals in the solution points
  Framework::DataSocketSink< CFreal > socket_solPntNormals;
  
  /// storage for the normals in the solution points
  Framework::DataSocketSink< CFreal > socket_flxPntNormals;
  
  /// storage for the cell volumes
  Framework::DataSocketSink< CFreal > socket_cellVolumes;
  
  /// storage of the geometrical jacobians in the states
  Framework::DataSocketSink<CFreal> socket_volumes;
  
  /// backup for epsilon
  CFreal m_epsBackUp;
  
  /// vector to store sol pnt values temporarily
  RealVector m_tempSolPntVec;
  
  /// vector to store sol pnt values temporarily
  RealVector m_tempSolPntVec2;
  
  /// boolean telling whether to print the LLAV output
  bool m_printLLAV;
  
  /// the gradients in the neighbouring cell
  std::vector< std::vector< std::vector< RealVector >* > > m_cellGradsAV;
  
  /// the corrected gradients in the flux points
  std::vector< std::vector< std::vector< RealVector* > > > m_cellGradFlxPntAV;
  
  /// average gradients in a flux point
  std::vector< RealVector* > m_avgGradAV;
  
  /// stores the flux jacobian for each side, in each sol pnt, for each variable, for each direction
  std::vector< std::vector< std::vector< std::vector< RealVector > > > > m_fluxJacobian;
  
  /// stores the Riemann flux jacobian for each side, in each face flx pnt, for each variable
  std::vector< std::vector< std::vector< RealVector > > > m_riemannFluxJacobian;
  
  /// convective Riemann flux
  std::vector < RealVector > m_flxPntRiemannFluxDiff;
  
  /// convective and diffusive Riemann flux (no LLAV)
  std::vector < RealVector > m_flxPntRiemannFluxDiffConv;
  
  /// perturbed convective Riemann Flux
  std::vector < RealVector > m_flxPntRiemannFluxPert;
  
  /// temporary storage for a flux
  RealVector m_tempFlux;
  
  /// stores the flux jacobian to the gradients for each side, in each sol pnt, for each variable, for each gradient direction for each flux direction
  std::vector< std::vector< std::vector< std::vector< std::vector< RealVector > > > > > m_gradientFluxJacobian;
  
  /// stores the gradient variables jacobian to the states for each side, in each sol pnt, for each depending state variable, for each grad vars variable
  std::vector< std::vector< std::vector< RealVector > > > m_gradVarsToStateJacobian;
  
  /// stores the gradient jacobian to the states for each side, in each sol pnt, for each depending side for each depending sol pnt, for each gradient direction
  std::vector< std::vector< std::vector< std::vector< RealVector > > > > m_gradientStateJacobian;
  
  /// stores the Riemann flux jacobian to the gradients for each face flx pnt, for each variable, for each gradient direction
  std::vector< std::vector< std::vector< RealVector > > > m_riemannFluxGradJacobian;
  
  /// Continuous diffusive flux at the solution points backup for both neighbor cells
  std::vector< std::vector< std::vector< RealVector > > > m_contFlxBackupDiff;
  
  /// variable to gradient variable dependency
  std::vector< std::vector< CFuint > > m_varToGradVarDep;
  
  /// number of gradient variable that depend on each variable
  std::vector< CFuint > m_nbrVarToGradVarDep;
  
  /// vector of booleans telling whether the current solution point needs to be added
  std::vector< bool > m_needToAddSolPnt;
  
  /// current faceID
  CFuint m_faceID;
  
  /// temporary state or state vector
  RealVector m_temp;
  
  /// temporary state or state vector for the other side
  RealVector m_tempOther;

  /// temporary state or state vector
  RealVector m_temp2;
  
  /// temporary state or state vector for the other side
  RealVector m_tempOther2;
  
  /// stores the AV jacobian to the state for each side, in each sol pnt, for each variable
  std::vector< std::vector< RealVector > > m_epsJacobian;
  
  /// stores the LLAV jacobian to the states for each side, in each sol pnt, for each variable to each side, to each sol pnt, to each flux direction
  std::vector< std::vector< std::vector< std::vector< std::vector< std::vector< RealVector > > > > > > m_llavFluxJacobian;
  
  /// stores the LLAV Riemann flux jacobian for each side, in each sol pnt, for each variable, for each flux point
  std::vector< std::vector< std::vector< std::vector< RealVector > > > > m_llavRiemannFluxJacobian;
  
  /// damping coefficient
  CFreal m_dampCoeffDiff;
  
  
  private:

  /// Physical data temporary vector
  RealVector m_pData;

}; // class ConvDiffLLAVJacobFluxReconstruction

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_ConvDiffLLAVJacobFluxReconstruction_hh
