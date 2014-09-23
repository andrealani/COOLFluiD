#ifndef COOLFluiD_Numerics_FluctSplit_STU_HOCRD_SplitStrategy_hh
#define COOLFluiD_Numerics_FluctSplit_STU_HOCRD_SplitStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/CFMat.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "FluctSplit/FluctuationSplitStrategy.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {
      class SpaceTime_Splitter;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a fluctuation splitting strategy of higher order
 * elements which residual is computed with a contour integral.
 *
 * @author Nadege Villedieu
 * @author Tiago Quintino
 *
 */
class STU_HOCRD_SplitStrategy : public FluctuationSplitStrategy {

public: // methods

  /**
   * Constructor.
   */
  STU_HOCRD_SplitStrategy(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~STU_HOCRD_SplitStrategy();

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
   * Compute the integral of the fluxes
   */
  void computeHOSTFluxIntegral(std::vector<Framework::State*>& states, std::vector<RealVector>& m_phi);
  /**
   * Compute the integral of the fluxes
   */
  void computeHOSTFluxIntegralP1(std::vector<Framework::State*>& states, std::vector<RealVector>& m_phi);

 /**
   * Compute the fluctuation for the first time iteration
   * @param residual the residual for each variable to distribute in each state
   */
  void docomputeFirstFluctuation(std::vector<RealVector>& residual);
  
 /**
   * Compute the fluctuation for the first time iteration
   * @param residual the residual for each variable to distribute in each state
   */
  void docomputeFirstFluctuationP1(std::vector<RealVector>& residual);


 /**
   * Compute the fluctuation for the first time iteration
   * @param residual the residual for each variable to distribute in each state
   */
  void docomputeFluctuation(std::vector<RealVector>& residual);

  void computeHOtime(CFuint i1, CFuint i2, CFuint i3);
  void computeHOtimeP1(CFuint i1, CFuint i2, CFuint i3);

  void setCurrentCell();
protected: // data
  /// Volume of the current cell
  CFreal CellVolume;

  /// The socket to use in this strategy for the update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;
  
  /// The socket to use in this strategy for the intermediate states
 Framework::DataSocketSink<Framework::State*> socket_interStates;

  /// temporary storage of the intermediate states
  std::vector<Framework::State*> m_interStates;

  /// the single splitter
  Common::SafePtr<SpaceTime_Splitter> m_splitter;

  /// Fluctuation for the intermediate layer of the sub-elements (advective part)
  std::vector<RealVector> m_phi_inter;

  /// Fluctruation for the present of the sub-elements (advective part)
  std::vector<RealVector> m_phi_present;

  /// states at quadrature points
  std::vector<Framework::State*> qdstates;

  /// extra vars at quadrature points
  std::vector<RealVector*> qdExtraVars;

  /// back up of update cell states
  std::vector<Framework::State*> m_statesBkp;

  /// temporary face normal
  RealVector facenormal;

  /// states that compose each sub face
  MathTools::CFMat<CFuint> subfacetable;

  /// quadrature points per face for advective fluctuation
  RealVector qd0a;
  RealVector qd1a;
  RealVector wqda;

  /// update variable set
  Common::SafePtr<Framework::ConvectiveVarSet> m_updateVar;

  /// fluxes on each sub face
  std::vector<RealVector> faceflux;

  /// coordinates of these states
  std::vector<Framework::Node*> qdnodes;

  /// direction of the faces
  MathTools::CFMat<CFreal> subelemfacedir;

  /// faces that compose each sub element
  MathTools::CFMat<CFuint> subelemtable;

  /// states that compose each sub-element
  MathTools::CFMat<CFuint> subelemtable_state;

  /// states that compose each sub element
  std::vector<Framework::State*> substates;

  /// Temporary states used to compute fluctuation from unsteady part
  std::vector<Framework::State*> temp_states;

  /// Fluctuation of the unsteadyness
  RealVector m_phi_time;

  /// quadrature points per face for time fluctuation
  RealVector qd0t;
  RealVector qd1t;
  RealVector qd2t;
  RealVector wqdt;

  /// temporary state
  RealVector m_u;

  /// Total fuctuation from the unsteadyness for the current sub-element
  RealVector m_flux_time;

  /// Total fuctuation from the advective part for the current sub-element
  RealVector m_flux_space;

  /// residuals for each sub element
  std::vector<RealVector> subresidual;

  /// The socket to use in this strategy for the past states
  Framework::DataSocketSink<Framework::State*> socket_pastStates;

  /// temporary storage of the past states
  std::vector<Framework::State*> m_pastStates;

  
  /// Fluctruation for the present of the sub-elements (advective part)
  std::vector<RealVector> m_phi_past;

  /// Storage for the residual (inter states contribution for the space residual)
  /// its size is nbcell and nbeqs*nbsubcell*nbstaInsubcell
  std::vector<RealVector> m_interResiduals;

    /// Storage for the residual (inter states contribution for the space residual)
  /// its size is nbcell and nbeqs*nbsubcell*nbstaInsubcell
  std::vector<RealVector> m_pastResiduals;

  /// Number of equations
  CFuint m_nbEqs;
  
  /// If this boolean is true the first time iteration is done as P1 in space
  /// otherwise it is P2 in space
  bool m_firstitP1;

}; // class STU_HOCRD_SplitStrategy

//////////////////////////////////////////////////////////////////////////////
      
    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_STU_HOCRD_SplitStrategy_hh
