#ifndef COOLFluiD_Framework_FVMCC_BaseBC_hh
#define COOLFluiD_Framework_FVMCC_BaseBC_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSink.hh"
#include "Framework/StateInterpolator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    
    class GeometricEntity;
    
//////////////////////////////////////////////////////////////////////////////

  /**
   * This functor class makes a comparison between two types
   *
   * @author Andrea Lani
   *
   */      
template <typename T1, typename T2>
class BiggerEqualVec {
public:
  bool operator()(const T1& v1, const T2& v2) {
    return (v1 >= v2);
  }
};
  
template <>
class BiggerEqualVec<RealVector, RealVector> {
public:
  bool operator()(const RealVector& v1, const RealVector& v2) 
  {
    cf_assert(v1.size() == v2.size());
    for (CFuint i = 0; i < v1.size(); ++i) {
      // each component of v1 has to be >= corresponding component of v2
      if (v1[i] < v2[i]) return false;
    }
    return true;
  }
};
      
//////////////////////////////////////////////////////////////////////////////
      
/**
 * This class represents the base class for any BC for cell center FV
 *
 * @author Andrea Lani
 * @modified by Alessandro Sanna
 *
 */
template <typename BASE>
class Framework_API FVMCC_BaseBC : public BASE {

public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  FVMCC_BaseBC(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~FVMCC_BaseBC();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();
  
  /**
   * Unsetup private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();
  
  /**
   * Set the preProcesses connectivity between faces belonging to different process
   *
   */
  virtual void preProcess();

 
  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Get the vector of flags specifying the variables that are 
   * constanly extrapolated (i.e. that use zero gradient)
   */
  std::vector<bool>* getZeroGradientsFlags() 
  {
    return &_zeroGradient;
  }

  /**
   * Apply boundary condition on the given face
   */
  virtual void setGhostState(Framework::GeometricEntity *const face) = 0;
  
  /**
   * Compute the flux in the current face
   */
  virtual void computeFlux(RealVector& result);

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
  /**
   * Set the flags telling if the ghost states have to be placed on the face itself
   */
  void setPutGhostsOnFace()
  { 
    _putGhostsOnFace = false;
    const std::vector<std::string>& trsNames = 
      this->getMethodData().getTRSsWithGhostsOnFace(); 
    for (CFuint i = 0; i < trsNames.size(); ++i) {
      if (trsNames[i] == this->getCurrentTRS()->getName()) {
	_putGhostsOnFace = true;
      }	
    }	
  }
  
  /// get the velocity IDs in the State (update variables)
  const std::vector<CFuint>& getStateVelocityIDs() const {return m_velocityIDs;}
  
protected: // functions

  /**
   * Execute on the current TRS
   */
  virtual void executeOnTrs();
    
  /**
   * Get the flags telling if the ghost states have to be placed on the face itself
   */
  bool getPutGhostsOnFace() const {return _putGhostsOnFace;}

  /**
   * Reposition the node if any of the given ghostState values is < 0
   */
  template <typename T1, typename T2>
  void repositionNode(const T1& innerValue, T1& ghostValue,
                      const T2& bValue, const T2& minValue);
  
  /**
   * This function makes a linear interpolation between the values in the
   * inner state and the ghost state ones
   */
  template <typename T1, typename T2, typename T3>
  void linearInterpolate(const T1& innerValue, const T2& wallValue, T3& ghostValue)
  {
    cf_assert(m_drXiXw > 0.);
    ghostValue = innerValue - (innerValue - wallValue)*(m_drXiXg/m_drXiXw);
  }
  
  /**
   * This function makes a linear interpolation between the values in the
   * inner state and the ghost state ones
   */
  template <typename T1, typename T2, typename T3>
  void linearInterpolateToB(const T1& innerValue, const T2& ghostValue, T3& wallValue)
  {
    cf_assert(m_drXiXw > 0.);
    wallValue = (ghostValue  - innerValue)/(m_drXiXg/m_drXiXw) + innerValue;
  }
  
   /**
    * Compute the original position of the ghost node
    */
  virtual void computeGhostPosition(Framework::GeometricEntity *const face);
  
  /// Get the @see StateInterpolator
  Common::SafePtr<StateInterpolator> getStateInterpolator() const 
  {
    return m_sInterpolator.getPtr();
  }
 
 private:
  
  /// state interpolator object
  Common::SelfRegistPtr<Framework::StateInterpolator> m_sInterpolator;
  
protected: // data
  
  /// put ghosts on face
  bool _putGhostsOnFace;
    
  /// handle to the face normals
  Framework::DataSocketSink< CFreal> socket_normals;

  /// handle to the face areas
  Framework::DataSocketSink< CFreal> socket_faceAreas;

  /// storage of the states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// storage of the ghost states
  Framework::DataSocketSink<Framework::State*> socket_gstates;

  /// storage of the nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;
  
  /// IDs corresponding to the velocity components
  std::vector<CFuint> m_velocityIDs;
  
  /// flag telling which variables to ignore
  std::vector<bool> m_computeVars;

  /// temporary coefficient  
  CFreal m_t;  
  
  /// temporary distance between the inner state and the wall  
  CFreal m_drXiXw;  
  
  /// temporary distance between the inner state and the repositioned ghost state  
  CFreal m_drXiXg;  
  
  /// factor  
  CFreal m_factor;  

  /// temporary internal node 
  RealVector* m_innerNode; 
  
  /// face mid point 
  RealVector m_faceMidPoint;
  
  /// vector (G-M) from the mid face M to the ghost point G
  RealVector m_XgXm;

  /// temporary node 
  RealVector m_tempNode; 
 
  /// temporary middle node 
  RealVector m_midNode; 
 
  /// temporary ghost node 
  RealVector m_tempGhostNode; 
 
  /// temporary ghost node backup 
  RealVector m_tempGhostNodeBkp; 

  /// temporary face normal 
  RealVector m_faceNormal; 
  
  /// map TRS name -> initial solution array that will be used as BC value
  Common::CFMap<std::string, RealVector*> m_initialSolutionMap;
  
  /// flag telling if a full BC loop over all the TRS faces is done
  bool m_fullLoop;
  
  /// vector of flags specifying the variables that are constanly extrapolated 
  /// (i.e. that use zero gradient)
  std::vector<bool> _zeroGradient;

  /// coefficient controlling the ghost node movement 
  CFreal m_coeff; 
  
  /// array specifying IDs of initial solution components that will be used as BC value
  std::vector<CFuint> m_initialSolutionIDs;
  
  /// name of state interpolator object
  std::string m_sInterpolatorStr;
  
}; // end of class FVMCC_BaseBC

//////////////////////////////////////////////////////////////////////////////

template <typename BASE>
template <typename T1, typename T2> 
void FVMCC_BaseBC<BASE>::repositionNode(const T1& innerValue, T1& ghostValue, 
					const T2& bValue, const T2& minValue) 
{ 
  using std::cout; 
  using std::endl; 
  using std::ofstream;

  BiggerEqualVec<T1, T2> testFun;
  ghostValue = 2.*bValue - innerValue; 
  m_factor = 1.0; 
  bool applyFix = false;  
  
  // CFuint count = 0; 
  // while (!(ghostValue >= minValue)) { 
  while (!testFun(ghostValue, minValue)) {
    applyFix = true;
    m_factor *= m_coeff; 
    // new position of the ghost node 
    m_tempNode = ((m_factor - 1.)*m_midNode  + m_tempGhostNode)/m_factor; 
    m_drXiXg = MathTools::MathFunctions::getDistance(*m_innerNode, m_tempNode); 
    // new temperature in the ghost state 
    linearInterpolate(innerValue, bValue, ghostValue); 
        
    // cout << "innerValue  = " << innerValue << endl; 
    // cout << "ghostValue  = " << ghostValue << endl; 
    // cout << "m_drXiXg = " << m_drXiXg << endl; 
    // cout << "m_drXiXw = " << m_drXiXw << endl; 
    // T ag((innerValue - ghostValue)/m_drXiXg); 
    // T aw((innerValue - m_wallTemp)/m_drXiXw); 
    // cout << "ag = " << ag << endl; 
    // cout << "aw = " << aw << endl; 
    // move the ghost to the new position 
    m_tempGhostNode = m_tempNode; 
  } 

  if (applyFix) { 
    // Check that the new ghost point is on the opposite side of the face
    // with respect to the inner node
    // Must be n*(G-M)>0  because n is outward and ghost node G should be outside 
    m_XgXm = m_tempGhostNode - m_midNode;
    const CFreal check = MathTools::MathFunctions::innerProd(m_faceNormal, m_XgXm);
    if (check < 0.) {
      cout << "### FVMCC_BaseBC::repositionNode() failed => n*(G-M) = " << check << " < 0" << endl;
      // string filename = "GhostWrong.txt" + StringOps::to_str(PE::GetPE().GetRank());
      // ofstream fout(filename);
      cout.precision(12); cout << "### inner T     =  " << innerValue << endl;
      cout.precision(12); cout << "### wall  T     =  " << bValue << endl;
      cout.precision(12); cout << "### face normal = (" << m_faceNormal << ")\n";
      cout.precision(12); cout << "### ghost node  = (" << m_tempGhostNode << ")\n";
      cout.precision(12); cout << "### projection  = (" << m_midNode << ")\n";
      cout.precision(12); cout << "### inner node  = (" << *m_innerNode << ")\n";
      //   fout.close();
    }
  } 

  cf_assert (ghostValue >= minValue);
} 
      
////////////////////////////////////////////////////////////////////////////// 

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "Framework/FVMCC_BaseBC.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_FVMCC_BaseBC_hh
