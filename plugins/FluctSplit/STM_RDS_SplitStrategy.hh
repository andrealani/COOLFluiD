#ifndef COOLFluiD_Numerics_FluctSplit_STM_RDS_SplitStrategy_hh
#define COOLFluiD_Numerics_FluctSplit_STM_RDS_SplitStrategy_hh

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
 * @author Thomas Wuilbaut
 *
 */
class STM_RDS_SplitStrategy : public FluctuationSplitStrategy {

public: // methods

  /**
   * Constructor.
   */
  STM_RDS_SplitStrategy(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~STM_RDS_SplitStrategy();

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
 /**
   * Sets the current cell and calls the computation of the
   * consistent state transformation.
   */
  void setCurrentCell();

private: // data
  
  /// The socket to use in this strategy for the update coefficients
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// The socket to use in this strategy for the past normals
  Framework::DataSocketSink< InwardNormalsData*> socket_pastNormals;

  /// The socket to use in this strategy for the past states
  Framework::DataSocketSink< Framework::State*> socket_pastStates;

  /// The socket for past volumes of the cells
  Framework::DataSocketSink<CFreal> socket_pastCellVolume;

  /// The socket for the speed of the cells
  Framework::DataSocketSink<RealVector> socket_cellSpeed;

  /// the single splitter
  Common::SafePtr<SpaceTime_Splitter> m_splitter;

  /// temporary storage of the past states
  std::vector<Framework::State*> _pastStates;

  ///Temporary vector
  RealVector m_tmpvec;

  /// Storage for the residual (past states contribution)
  std::vector<RealVector> _pastResiduals;

  /// Storage for the residual (past states contribution of the N scheme used in B scheme)
  std::vector<RealVector> _pastResiduals_order1;

  /// temporary vector with past source term residual
  std::vector<RealVector> temp_residual;

/// back up of update cell states
  std::vector<Framework::State*> m_statesBkp;
}; // class STM_RDS_SplitStrategy

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_STM_RDS_SplitStrategy_hh
