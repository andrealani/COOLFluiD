#ifndef COOLFluiD_Numerics_FluctSplit_STM_HOCRDS2T2_SplitStrategy_hh
#define COOLFluiD_Numerics_FluctSplit_STM_HOCRDS2T2_SplitStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctSplit/FluctuationSplitStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {
      class SpaceTime_Splitter;
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a fluctuation splitting strategy
 *
 * @author Nadege Villedieu
 *
 */
class STM_HOCRDS2T2_SplitStrategy : public FluctuationSplitStrategy {

public: // methods

  /**
   * Constructor.
   */
  STM_HOCRDS2T2_SplitStrategy(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~STM_HOCRDS2T2_SplitStrategy();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    FluctuationSplitStrategy::configure(args);
  }

  /**
   * Set up private data and data
   */
  virtual void setup();

  /**
   * Unsetup the private data and data of the aggregated classes
   * in this strategy after the  processing phase
   */
  virtual void unsetup();

  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Compute the fluctuation
   * @param residual the residual for each variable to distribute in each state
   */
  virtual void computeFluctuation(std::vector<RealVector>& residual);

  static void defineConfigOptions(Config::OptionList& options);

  protected: // methods

  /**
   * Compute the present contribution of integral of the Space-Time fluxes 
   */

  void computeSTFluxIntegral();

/**
   * Compute the past contribution of integral of the Space-Time fluxes
   */

  void computeHOFluctuation(std::vector<Framework::State*> & states);
 
 /**
   * Sets the current cell and calls the computation of the
   * consistent state transformation.
   */

 /**
   * Compute the fluctuation for the first time iteration
   * @param residual the residual for each variable to distribute in each state
   */
  void docomputeFirstFluctuation(std::vector<RealVector>& residual);

 /**
   * Compute the fluctuation for all the others time iteration
   * @param residual the residual for each variable to distribute in each state
   */
  void docomputeFluctuation(std::vector<RealVector>& residual);

  void setCurrentCell();

private: // data
  
  /// The socket to use in this strategy for the update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// The socket to use in this strategy for the past normals
  Framework::DataSocketSink<InwardNormalsData*> socket_pastNormals;

  /// The socket to use in this strategy for the past states
  Framework::DataSocketSink<Framework::State*> socket_pastStates;

  /// The socket to use in this strategy for the past states
  Framework::DataSocketSink<Framework::State*> socket_interStates;

  /// The socket for past volumes of the cells
  Framework::DataSocketSink<CFreal> socket_pastCellVolume;

  /// The socket for the speed of the cells
  Framework::DataSocketSink<RealVector> socket_cellSpeed;

  /// the single splitter
  Common::SafePtr<SpaceTime_Splitter> m_splitter;

  /// temporary storage of the past states
  std::vector<Framework::State*> _pastStates;

  /// temporary storage of the past states
  std::vector<Framework::State*> _interStates;

  /// temporary storage of the consistent past states
  std::vector<Framework::State*> * tPastStates;

  ///Temporary vector
  RealVector _null;

  /// Contribution from the past to the residual of the cell
  RealVector m_flux_past_time;

  /// Contribution from the present to the residual of the cell
  RealVector m_flux_time;

  /// Contribution from the past to the residual of the cell
  RealVector m_flux_past_space;

  /// Contribution from the present to the residual of the cell
  RealVector m_flux_space;

  /// Volume of the current cell
  CFreal CellVolume;

 /// interpolated states at quadrature points
  std::vector<Framework::State*> m_qdstates;

/// extra values at quadrature points
  std::vector<RealVector*> m_qdExtraVars;

  /// contour integrator
  Common::SafePtr<Framework::ContourIntegrator> m_contourIntegrator;

  /// array for the number of quadrature points in each cell
  /// @todo maybe this could be handled on th side of the integrator itself
  std::vector<CFuint> m_nbQPointsInCell;

  /// matrix storing the unit sized face normals
  std::vector<RealVector> m_unitFaceNormals;

  /// back up of update cell states
  std::vector<Framework::State*> m_statesBkp;

  ///Temporary vector of the contour integration
  RealVector m_phi;

  /// Storage for the residual (past states contribution)
  std::vector<RealVector> _pastResiduals;

  /// Storage for the residual (past states contribution of the N scheme used in B scheme)
  std::vector<RealVector> _pastResiduals_order1;

  /// temporary vector with past source term residual
  std::vector<RealVector> temp_residual;

/// fluctuation on each sub element
  std::vector<RealVector*> m_phisubT;
  
  /// direction of the faces
  MathTools::CFMat<CFreal> subelemfacedir;

  /// faces that compose each sub element
  MathTools::CFMat<CFuint> subelemtable;

  /// states that compose each sub face
  MathTools::CFMat<CFuint> subfacetable;

  /// states that compose each sub element
  std::vector<Framework::State*> substates;

  /// residuals for each sub element
  std::vector<RealVector> subresidual;

  /// fluxes on each sub face
  std::vector<RealVector> faceflux;

  /// quadrature points per face
  RealVector qd0;
  RealVector qd1;
  RealVector wqd;

  /// states at quadrature points
  std::vector<Framework::State*> qdstates;
  
  /// extra vars at quadrature points
  std::vector<RealVector*> qdExtraVars;
  
  /// coordinates of these states
  std::vector<Framework::Node*> qdnodes;

  /// temporary face normal
  RealVector facenormal;


  /// update variable set
  Common::SafePtr<Framework::ConvectiveVarSet> m_updateVar;

/// Coeficients used to compute the mass matrix
  MathTools::CFMat<CFreal> kappa1;
  MathTools::CFMat<CFreal> kappa2;
  MathTools::CFMat<CFreal> kappa3;
  MathTools::CFMat<CFreal> kappa4;
  MathTools::CFMat<CFreal> kappab1;
  MathTools::CFMat<CFreal> kappab2;
  MathTools::CFMat<CFreal> kappab3;
  MathTools::CFMat<CFreal> kappab4;

  /// Temporary vectors
  RealVector temp_vect0;
  RealVector temp_vect1;
  RealVector temp_vect2;
  RealVector temp_vect3;
  RealVector temp_vect4;
  RealVector temp_vect5;

  /// This coeffisients are used to define the time discretization :
  /// (1+m_xi)*states - (1+2m_xi)*interSate + m_xi*pastState = 
  ///                                dt*(m_theta*res + (1-m_theta-m_phi)*res_inter - m_phi*res_past)
  CFreal m_theta;
  CFreal m_xi;
  CFreal m_chi;

}; // class STM_HOCRDS2T2_SplitStrategy

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_STM_HOCRDS2T2_SplitStrategy_hh
