// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_LLAVFluxReconstruction_hh
#define COOLFluiD_FluxReconstructionMethod_LLAVFluxReconstruction_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

#include "FluxReconstructionMethod/DiffRHSFluxReconstruction.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/// Command to add Localized Laplacian Artificial Viscosity near discontinuities
/// @author Ray Vandenhoeck
    
class LLAVFluxReconstruction : public DiffRHSFluxReconstruction {

public: // functions

  /// Constructor
  explicit LLAVFluxReconstruction(const std::string& name);

  /// Destructor
  virtual ~LLAVFluxReconstruction() {}

  /// Execute processing actions
  void execute();
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unsetup private data
   */
  virtual void unsetup();
    
protected: //functions
  
  /**
   * compute the wave speed updates for this face
   * @pre reconstructFluxPntsStates(), reconstructFaceAvgState(),
   *      setFaceTermData() and set the geometrical data of the face
   */
  virtual void computeWaveSpeedUpdates(std::vector< CFreal >& waveSpeedUpd);

  /// compute the interface flux
  void computeInterfaceFlxCorrection();
  
  /// compute the divergence of the discontinuous flux (-divFD+divhFD)
  void computeDivDiscontFlx(std::vector< RealVector >& residuals);
  
  /**
   * Set the data for the current face necessary to calculate FI
   */
  virtual void setFaceData(CFuint faceID);
  
  /**
   * Set the data for the current cell necessary to calculate the residual update
   */
  virtual void setCellData();
  
  /**
   * Compute the projected states on order P-1
   */
  void computeProjStates(std::vector< RealVector >& projStates);
  
  /**
   * Compute the artificial viscosity
   */
  void computeEpsilon();
  
  /**
   * Compute the reference artificial viscosity
   */
  void computeEpsilon0();
  
  /**
   * Compute the smoothness indicator
   */
  void computeSmoothness();
  
  /**
   * Store the computed artificial viscosity
   */
  void storeEpsilon();
  
  /**
   * Add extra viscosity to regions where positivity is almost violated, using the minValue option comparing to the state of the monitoredVar option
   */
  void addPositivityPreservation();
  
  /**
   * command to compute the peclet number based on user input
   */
  virtual CFreal computePeclet(); 
  
  /**
   * command to compute the coefficient to be multiplied with epsilon for the wavespeedupdates
   */
  virtual CFreal computeViscCoef(RealVector* state);

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
  std::vector< CFreal > m_solEpsilons;
  
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
  
  /// peclet number
  CFreal m_peclet;
  
  /// average artificial viscosities in the nodes
  RealVector m_nodeEpsilons;
  
  /// number of neighbors for each node
  RealVector m_nbNodeNeighbors;
  
  /// average artificial viscosities in the elements
  RealVector m_cellEpsilons;
  
  /// vector containing pointers to the nodes in a cell
  std::vector< Framework::Node*  >* m_cellNodes;
  
  /// vector containing pointers to the nodes in a face
  std::vector< Framework::Node*  >* m_faceNodes;
  
  /// number of corner nodes for current element type
  CFuint m_nbrCornerNodes;
  
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
  
  /// builder of cells
  Common::SafePtr< Framework::GeometricEntityPool<CellToFaceGEBuilder> > m_cellBuilder;
  
  /// pointer to booleans telling whether a face is on the boundary
  Common::SafePtr< std::vector< bool > > m_isFaceOnBoundaryCell;

  /// pointer to neighbouring cell side vector
  Common::SafePtr< std::vector< CFuint > > m_nghbrCellSideCell;

  /// pointer to current cell side vector
  Common::SafePtr< std::vector< CFuint > > m_currCellSideCell;

  /// pointer to orientation vector
  Common::SafePtr< std::vector< CFuint > > m_faceOrientsCell;

  /// pointer to BC index vector
  Common::SafePtr< std::vector< CFuint > > m_faceBCIdxCell;
  
  /// variable for faces
  const std::vector< Framework::GeometricEntity* >* m_faces;
  
  /// boundary condition state computers
  Common::SafePtr< std::vector< Common::SafePtr< BCStateComputer > > > m_bcStateComputers;
  
  /// the ghost gradients in the flux points
  std::vector< std::vector< RealVector* > > m_flxPntGhostGrads;
  
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
  
  /// bool telling whether the Jacobian is being computed
  bool m_jacob;
  
  /// bool telling whether to add viscosity for positivity preservation
  bool m_addPosPrev;
  
  /// boolean telling whether to add the contribution of the artificial flux to the update coefficients
  bool m_addUpdCoeff;
  
  /// value at which point positivity preservation is added
  CFreal m_minValue;
  
  /// index of the monitored variable for positivity preservation
  CFuint m_monitoredVar;
  
  /// number of times positivity preservation activated
  CFuint m_nbPosPrev;
  
  /// number of times positivity preservation activated over all CPUs
  CFuint m_nbPosPrevGlobal;
  
  /// maximum factor applied to viscosity for positivity preservation
  CFreal m_viscFactor;
  
  /// subcell resolution
  CFreal m_subcellRes;
  
  private:

  /// Physical data temporary vector
  RealVector m_pData;
  
}; // class Solve

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_LLAVFluxReconstruction_hh

