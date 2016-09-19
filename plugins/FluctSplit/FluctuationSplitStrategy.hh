#ifndef COOLFluiD_Numerics_FluctSplit_FluctuationSplitStrategy_hh
#define COOLFluiD_Numerics_FluctSplit_FluctuationSplitStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/MethodStrategy.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "Framework/BaseDataSocketSink.hh"
#include "FluctSplit/FluctuationSplitData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework { class GeometricEntity; }



      namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represent a fluctuation splitting strategy
/// @author Andrea Lani
/// @author Tiago Quintino
class FluctSplit_API FluctuationSplitStrategy : public Framework::MethodStrategy<FluctuationSplitData> {

public: // typedefs

  typedef Framework::BaseMethodStrategyProvider<FluctuationSplitData, FluctuationSplitStrategy> PROVIDER;
  typedef Framework::VarSetMatrixTransformer MatrixTransformer;
  typedef Framework::VarSetTransformer       VectorTransformer;

public: // methods

  /// Constructor
  FluctuationSplitStrategy(const std::string& name);

  /// Destructor
  virtual ~FluctuationSplitStrategy();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure ( Config::ConfigArgs& args );

  /// Sets up the method data
  virtual void setup();

  /// Unsets the method data
  virtual void unsetup();

  /// function called at each prepare phase before the computation of the RHS
  virtual void prepare ();

  /// Returns the DataSocket's that this numerical strategy needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Computes the fluctuation
  virtual void computeFluctuation(std::vector<RealVector>& residual) = 0;

  /// Gets the Class name
  static std::string getClassName() { return "FluctuationSplitStrategy"; }
  
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
protected: // methods
  
  /// Compute the physical data in the states
  /// This function is inlinable and should remain here
  void computeStatesData(CFuint nbPoints, 
			 Common::SafePtr<Framework::ConvectiveVarSet> vs,
			 const std::vector<Framework::State*>& qstates, 
			 std::vector<RealVector>& pdata,
			 std::vector<RealVector*>& extraVars)
  {
    for (CFuint i = 0; i < nbPoints; ++i) {
      vs->setExtraPhysicalVars(extraVars[i]);
      vs->computePhysicalData(*qstates[i], pdata[i]); 
    }
  }
  
  /// Compute the transformed states
  /// @pre assumes m_cell has been set to the current geometric entity
  /// @param states the states to be transformed to the consistent states
  /// @return the transformed consistent states
  std::vector<Framework::State*>*
  computeConsistentStates (std::vector<Framework::State*> *const states);

protected: // data

  /// The sockets to use in this strategy for the normals
  Framework::DataSocketSink< InwardNormalsData*> socket_normals;

  /// The sockets for the flags telling if states are on the boundary
  Framework::DataSocketSink<bool> socket_isBState;

  /// The sockets for the update coefficients
  Framework::DataSocketSink<CFreal> socket_volumes;

  /// handle for the InnerCells trs
  Common::SafePtr<Framework::TopologicalRegionSet> _cells;

  /// transformed linear states in current cell
  std::vector<Framework::State*> * m_linearStates;

  /// builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> _stdTrsGeoBuilder;
  
  /// interpolated physical data at quadrature points
  std::vector<RealVector> m_pdata;
  
}; // class FluctuationSplitStrategy

//////////////////////////////////////////////////////////////////////////////

inline std::vector<Framework::State*> *
FluctuationSplitStrategy::computeConsistentStates(std::vector<Framework::State*> *const states)
{  
  FluctuationSplitData& mdata = getMethodData();
  // transform the update states to linear varset
  m_linearStates = mdata.getUpdateToLinearVecTrans()->transform(states);
  // computes an average state in which jacobians will be linearized
  mdata.getLinearizer()->linearize(*m_linearStates);
  // transformation from linear to consistent variables
  return mdata.getLinearToDistribMatTrans()->transformFromRef(m_linearStates);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_FluctuationSplitStrategy_hh
