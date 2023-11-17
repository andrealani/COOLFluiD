#ifndef COOLFluiD_Framework_NodalStatesExtrapolator_hh
#define COOLFluiD_Framework_NodalStatesExtrapolator_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/OwnedObject.hh"
#include "Common/SafePtr.hh"
#include "Common/CFMap.hh"
#include "Environment/ConcreteProvider.hh"
#include "Framework/Storage.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/MethodStrategy.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/StateInterpolator.hh"
#include "MathTools/RealVector.hh"
#include "MathTools/RealMatrix.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework
  {
    class Node;
    class State;
    class TopologicalRegionSet;

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a nodal states extrapolator object
 *
 * @author Andrea Lani
 *
 */
template <typename DATA>
class NodalStatesExtrapolator : public Framework::MethodStrategy<DATA> {
public:

  typedef Framework::BaseMethodStrategyProvider<DATA,NodalStatesExtrapolator<DATA> > PROVIDER;
  
  typedef Common::CFMap<Framework::TopologicalRegionSet*, RealVector*> MapTrs2NodalValues;
  typedef Common::CFMap<Framework::TopologicalRegionSet*, Common::CFMap<CFuint,CFuint>*> MapTrsNodeIDs;

  /**
   * Constructor
   */
  NodalStatesExtrapolator(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NodalStatesExtrapolator();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Set up private data needed by the computation
   */
  virtual void setup(); 
  
  /**
   * Unsetup private data needed by the computation
   */
  virtual void unsetup();

  /**
   * Extrapolate the solution in all mesh nodes
   */
  virtual void extrapolateInAllNodes() = 0;

  /**
   * Extrapolate the solution in the given nodes
   */
  virtual void extrapolateInNodes
  (const std::vector<Framework::Node*>& nodes) = 0;

  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result;

    result.push_back(&socket_nodes);
    result.push_back(&socket_nstates);
    result.push_back(&socket_trsID);
    result.push_back(&socket_states);
    result.push_back(&socket_gstates);

    return result;
  }

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "NodalStatesExtrapolator";
  }
  
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
  /**
   * Gets list of state neighbors for a given node
   */
  const std::vector<Framework::State*>& getNodalStateNeighbors(CFuint nodeID) const
  {
    return _neighborStates[nodeID];
  }

  /**
   * Get the mapping TopologicalRegionSet - boundary nodal values
   */
  Common::SafePtr<std::vector<MapTrs2NodalValues*> > getMapTrs2NodalValues()
  {
    return &_mapTrs2NodalValues;
  }

  /**
   * Get the mapping TopologicalRegionSet - local nodeIDs to local TRS node IDs
   */
  Common::SafePtr<MapTrsNodeIDs> getMapTrs2NodeIDs()
  {
    return &_mapNodeID2TrsNodeID;
  }
  
  /**
   * Get nodal value corresponding to the given TRS, variableID and
   * local (in the current mesh partition) node ID
   */
  CFreal getNodalValue(Common::SafePtr<Framework::TopologicalRegionSet> trs,
		       CFuint variableID, CFuint localNodeID)
  {
    cf_assert(trs.isNotNull());
    const RealVector& valuesArray = *_mapTrs2NodalValues[variableID]->find(&*trs);
    return valuesArray[_mapNodeID2TrsNodeID.find(&*trs)->find(localNodeID)];
  }
  
  
  /// Add an array of lookup tables associated to a particular TRS 
  /// @pre to be filled in during setup() by specific BC Commands
  void addLookupState(const std::string& trsName,
		      Common::SafePtr<StateInterpolator> interp)
  {
    m_trsID2LookupState.insert(_mapTrsNameToID.find(trsName), interp);
  }
  
  /// return the number of iterations to run adiabatic
  bool runAdiabatic() const 
  {return (Framework::SubSystemStatusStack::getActive()->getNbIter() < m_nbIterAdiabatic);} 

  /**
   * Extrapolate some variables from a given file in time
   */
  virtual void extrapolateVarsFromFileInTime();
  
protected: // helper function
  
  /**
   * Local nested class for grouping surface data
   */
  class SurfaceData {
  public:
    RealMatrix xyz;
    RealVector Tw;
  };
  
  /**
   * Local nested class for grouping data related to the closet point
   */
  class ClosestPointData {
  public:
    std::vector<CFint> surfaceIDs;
    std::vector<CFint> pointsIDs;
    RealVector r;
    
    void reset() 
    {
      surfaceIDs.assign(surfaceIDs.size(), -1); 
      pointsIDs.assign(pointsIDs.size(), -1);
      r = MathTools::MathConsts::CFrealMax();
    }
    
    void regressionFromTo(CFuint start, CFuint end)
    {
      cf_assert(end < surfaceIDs.size());
      cf_assert(start < surfaceIDs.size());
      surfaceIDs[end] = surfaceIDs[start]; 

      cf_assert(end < pointsIDs.size());
      cf_assert(start < pointsIDs.size());
      pointsIDs[end] = pointsIDs[start];

      cf_assert(end < r.size());
      cf_assert(start < r.size());
      r[end] = r[start];
    }
  };
    
  /**
   * Extrapolate some variables from a given file
   */
  void extrapolateVarsFromFile(const std::vector<SurfaceData*>& surfaces);
  
  /**
   * Allocate mapping data needed for interpolation 
   */
  void allocateMappingData();
  
  /**
   * Read the surface data
   */
  void readSurfaceData(std::vector<SurfaceData*>& surfaces,
		       const std::string& fileName);
  
  /**
   * Read line data at z=0 from given surface
   */
  SurfaceData* extractLineData(SurfaceData* surface);
  
  /**
   * Add the boundary neighbors to the list of the state neighbors
   * for each node
   */
  virtual void addBoundaryNeighbors
  (std::vector< std::vector<CFuint> >& bFacesPerNode,
   Common::CFMap<CFuint, CFint>& mapFaceTrsID);

protected:

  /// socket for nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// socket for nodal States
  Framework::DataSocketSink< RealVector> socket_nstates;

  /// socket for TRS ID's (they are set according to _trsPriorityList and _orderedTrsList)
  Framework::DataSocketSink< CFint> socket_trsID;

  /// socket for states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// socket for ghost states
  Framework::DataSocketSink< Framework::State*> socket_gstates;
  
  /// current node ID
  CFuint _currNodeID;
  
  /// array of flags telling which interpolated or previously computed nodal values are available
  std::vector<bool> _nodalValuesIDFlags;

  /// map each TopologicalRegionSet with corresponding boundary nodal values
  std::vector<MapTrs2NodalValues*> _mapTrs2NodalValues;

  // reverse mapping between a local node ID and local node ID inside a specific TRS
  MapTrsNodeIDs _mapNodeID2TrsNodeID;
  
  /// map the boundary TRS names with a unique TRS ID
  Common::CFMap<std::string, CFint> _mapTrsNameToID;
  
  /// look up table for nodal values
  Common::CFMap<CFuint, Common::SafePtr<StateInterpolator> > m_trsID2LookupState;
  
  /// list of TRSs which has been reordered by priority
  std::vector<Common::SafePtr<Framework::TopologicalRegionSet> > _orderedTrsList;

  /// storage of the constant lcoefficients
  std::vector<std::vector<Framework::State*> > _neighborStates;
  
  /// storage of surface data
  std::vector<std::vector<SurfaceData*> > m_allSurfaces;
  
  /// storage of surface data at given time
  std::vector<SurfaceData*> m_surfaceAtTime;
  
  /// ID corresponding to the z coordinate (z=0) for which plane is extracted
  CFint m_extractCoordZID;
  
  /// names of the TRS defining the priority list
  std::vector<std::string> _trsPriorityList;
  
  // name of the TRSs on which values must be prescribed
  std::vector<std::string> _trsName;
  
  /// name of the file where the boundary distribution is provided
  std::vector<std::string> m_fileNameTw;
  
  /// time corresponding to the file where the boundary distribution is provided
  std::vector<CFreal> m_fileNameTime;
  
  /// rotation angle
  CFreal m_angle;
  
  /// ID  of the spatial coordinates lying in the rotation plane
  std::vector<CFuint> m_xvec;
  
  /// ID of the temperature
  CFuint m_tempID;
  
  /// number of closest points for interpolation stencil
  CFuint m_nbClosestPoints;
  
  /// IDs corresponding to the x,y coordinate (z=0) for which plane is extracted
  std::vector<CFint> m_extractCoordXYID;

  /// number of iterations to run adiabatic
  CFuint m_nbIterAdiabatic;
  
}; // end of class NodalStatesExtrapolator
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "Framework/NodalStatesExtrapolator.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_NodalStatesExtrapolator_hh
