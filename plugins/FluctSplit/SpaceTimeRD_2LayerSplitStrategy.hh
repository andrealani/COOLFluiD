#ifndef COOLFluiD_Numerics_FluctSplit_SpaceTimeRD_2LayerSplitStrategy_hh
#define COOLFluiD_Numerics_FluctSplit_SpaceTimeRD_2LayerSplitStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluctuationSplitStrategy.hh"

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
class SpaceTimeRD_2LayerSplitStrategy : public FluctuationSplitStrategy {

public: // methods

  /**
   * Constructor.
   */
  SpaceTimeRD_2LayerSplitStrategy(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~SpaceTimeRD_2LayerSplitStrategy();

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
  
private:
  
  /**
   * Compute the fluctuation and the update coefficient
   * @param residual the residual for each variable to distribute in each state
   */
  virtual void doComputeFluctAndUpdateCoeff(std::vector<RealVector>& residual);
  
  /**
   * Compute the fluctuation and the update coefficient
   * @param residual the residual for each variable to distribute in each state
   */
  virtual void doComputeFluct(std::vector<RealVector>& residual);
    
private: // data
  
  /// The socket to use in this strategy for the update coefficient
  Framework::DataSocketSink<CFreal> socket_updateCoeff;

  /// The socket to use in this strategy for the past normals
  Framework::DataSocketSink< InwardNormalsData*> socket_pastNormals;

  /// The socket to use in this strategy for the past normals (intermediate)
  Framework::DataSocketSink< InwardNormalsData*> socket_interNormals;

  /// The socket to use in this strategy for the past states
  Framework::DataSocketSink< Framework::State*> socket_pastStates;

  /// The socket to use in this strategy for the past states
  Framework::DataSocketSink< Framework::State*> socket_interStates;

  /// The socket for past volumes of the cells
  Framework::DataSocketSink<CFreal> socket_pastCellVolume;

  /// The socket for the speed of the cells
  Framework::DataSocketSink<RealVector> socket_cellSpeed;

  /// The socket for update coef for intermediate states
  Framework::DataSocketSink<CFreal> socket_interUpdateCoeff;

  /// the single splitter
  Common::SafePtr<SpaceTime_Splitter> m_splitter;

  /// temporary storage of the past states
  std::vector<Framework::State*> _pastStates;

  /// temporary storage of the intermediate states
  std::vector<Framework::State*> _interStates;

  ///Temporary vector
  RealVector _null;

  /// Storage for the residual (past states contribution)
  std::vector<RealVector> _pastResiduals;

  /// Storage for the residual (intermediate states)
  std::vector<RealVector> _interResiduals;

  /// Temporary storage for the intermediate Volume
  CFreal _interVolume;

}; // class SpaceTimeRD_2LayerSplitStrategy

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_SpaceTimeRD_2LayerSplitStrategy_hh
