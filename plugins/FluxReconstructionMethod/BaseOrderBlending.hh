// Copyright (C) 2016 KU Leuven, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#ifndef COOLFluiD_FluxReconstructionMethod_BaseFluxFiltering_hh
#define COOLFluiD_FluxReconstructionMethod_BaseFluxFiltering_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"

#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command to compute the blending coefficient for the order blending
 *
 * @author Rayan Dhib
 *
 */
class BaseOrderBlending : public FluxReconstructionSolverCom {
public:

  /**
   * Constructor.
   */
  explicit BaseOrderBlending(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~BaseOrderBlending();

  /**
   * Setup private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unsetup private data
   */
  virtual void unsetup();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  //virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > >
  //  needsSockets();
    
  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSource > >
    providesSockets();

  /**
   * Execute Processing actions
   */
  void execute();

protected: // functions

  /**
  * Compute the smoothness indicator
  */
  virtual void computeSmoothness();

  /**
  * Compute the smoothness indicator using modal energy method
  */
  virtual void computeSmoothnessModal();

  /**
  * Compute the smoothness indicator using projection method
  */
  virtual void computeSmoothnessProjection();

  /**
  * Compute the projected states on order P-1
  */
  void computeProjStates(std::vector< RealVector >& projStates);

  /**
  * [DEPRECATED] Compute threshold for blending coefficient (sigmoid method)
  * T = a * 10^(-c * (N+1)^0.25)
  */
  CFreal computeThreshold_Legacy(CFuint N, CFreal a, CFreal c);

  /**
  * Compute reference threshold for log-scale method
  * s0 = -S0 * log10(N+1)
  */
  CFreal computeReferenceThreshold() const;

  /**
  * Compute blending coefficient from smoothness indicator
  * If m_useLogScale=true: sinusoidal ramp in [s0-kappa, s0+kappa]
  * If m_useLogScale=false: [DEPRECATED] sigmoid function
  */
  CFreal computeBlendingCoefficient(CFreal smoothness);

  /**
  * Apply alpha limits to blending coefficient
  */
  CFreal applyAlphaLimits(CFreal alpha);

  /**
   * Computes vector of maximum modal order
   */
  RealVector getmaxModalOrder(const CFGeoShape::Type elemShape, const CFuint m_order);
  

protected: // data
  
  /// socket for output of the blending coefficients
  Framework::DataSocketSource< CFreal > socket_alpha;

  /// socket for previous blending coefficients
  Framework::DataSocketSource< CFreal > socket_prevAlpha;

  /// storage for the smoothness
  Framework::DataSocketSource<CFreal> socket_smoothness;

  /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> > m_cellBuilder;

    /// builder of cells
  Common::SafePtr<Framework::GeometricEntityPool<CellToFaceGEBuilder> > m_geoBuilder;

  /// variable for cell
  Framework::GeometricEntity* m_cell;

  /// vector containing pointers to the states in a cell
  std::vector< Framework::State* >* m_cellStates;

  /// current element index
  CFuint m_elemIdx;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

  /// index of element type
  CFuint m_iElemType;

  /// dimensionality of the physical model
  CFuint m_dim;
  
  /// number of solution points for current element type
  CFuint m_nbrSolPnts;
  
  /// maximum number of flux points in a cell
  CFuint m_maxNbrFlxPnts;
  
  /// order of FR method
  CFuint m_order;
  
  /// coefs to extrapolate the states to the flx pnts
  Common::SafePtr< std::vector< std::vector< CFreal > > > m_solPolyValsAtFlxPnts;
  
  /// extrapolated states in the flux points of the cell
  std::vector< RealVector > m_cellStatesFlxPnt;
  
  /// number of times filtering was applied
  CFuint m_nbLimits;
  
  /// total number of times a cell was filtered for all processors
  CFuint m_totalNbLimits;
  
  /// showrate of PP info
  CFuint m_showrate;

  /// iteration at which the filter strength will be frozen
  CFuint m_freezeFilterIter;

  /// IDs of face neighboring elements
  std::vector< std::vector<CFuint> > m_NeighborIDs;
  
  /// Element-wise filter strength
  CFreal m_filterStrength;

  /// vector of temporary cell states
  std::vector< RealVector > m_tempStates;

  /// Vandermonde Matrix
  RealMatrix  m_vdm;
  
  /// Inverse of the Vandermonde
  RealMatrix m_vdmInv;

  /// Vector of maximum modal order
  RealVector m_maxModalOrder;

  /// Nb of solution points for P-1
  CFuint m_nbrSolPntsMinOne;

  /// Nb of solution points for P-2
  CFuint m_nbrSolPntsMinTwo; 

  /// reference smoothness
  CFreal m_s0;

  /// smoothness
  CFreal m_s;

  /// maximum smoothness in domain
  CFreal m_Smax;
  
  /// controlling parameter kappa
  CFreal m_kappa;

  /// index of the monitored variable for positivity preservation
  CFuint m_monitoredVar;

  /// states projected on P-1
  std::vector< RealVector > m_statesPMinOne;

  /// transformation matrices to order P-1
  RealMatrix m_transformationMatrix;

  /// vector to store sol pnt values temporarily
  RealVector m_tempSolPntVec;

  /// vector to store sol pnt values temporarily
  RealVector m_tempSolPntVec2;

  // Order blending specific parameters
  /// Number of smoothing sweeps
  CFuint m_nbSweeps;

  /// Threshold function parameter 'a'
  CFreal m_thresholdA;

  /// Threshold function parameter 'c'
  CFreal m_thresholdC;

  /// Sigmoid function parameter 's'
  CFreal m_sigmoidS;

  /// Minimum blending coefficient
  CFreal m_alphaMin;

  /// Maximum blending coefficient
  CFreal m_alphaMax;

  /// Weight factor for neighbor influence
  CFreal m_neighborWeight;

  /// Weight factor for neighbor influence during sweeps
  CFreal m_sweepWeight;

  /// Damping factor for sweep iterations
  CFreal m_sweepDamping;

  /// Use log-scale smoothness indicator (true, recommended) or legacy sigmoid (false)
  bool m_useLogScale;

  /// Method for computing smoothness ("Modal" or "Projection")
  std::string m_smoothnessMethod;

  /// Expression for monitored variable in modal method
  std::string m_modalMonitoredExpression;

}; // class BaseOrderBlending

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_BaseFluxFiltering_hh