#ifndef COOLFluiD_Numerics_ExplicitFilters_InterpolationStencilComputer_hh
#define COOLFluiD_Numerics_ExplicitFilters_InterpolationStencilComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "ExplicitFilters/StencilComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace ExplicitFilters {
            
      class Triangle;

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class provides an abstract interface for functors
   * computing the stencil for a specific kind of polynomial
   * reconstruction in ExplicitFilters
   *
   * @author Willem Deconinck
   */
class InterpolationStencilComputer : public StencilComputer {

public:

  /**
   * Returns the DataSockets that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> >
  providesSockets();

  /**
   * Defines the config options of this class
   * @param   options   config options of this class
   */
  static void defineConfigOptions(Config::OptionList& options);

  /** 
   * Configure the object
   * @param    args   configuration arguments
   */
  virtual void configure(Config::ConfigArgs& args);

  /**
   * Constructor
   */
   InterpolationStencilComputer(const std::string& name);

  /**
   * Destructor
   */
  virtual ~InterpolationStencilComputer();

  /**
   * Setup
   */
  virtual void setup();

  /**
   * Post-process a specific stencil
   */
  virtual void postProcessStencil(const CFuint& centreStateID);

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "InterpolationStencilComputer";
  }

  CFuint CFround(const CFreal& x) {
    return floor(x+0.5);
  }

protected: //data

  /// storage for the stencil via triangles
  Framework::DataSocketSource < std::vector<Triangle> > socket_triangles;
  
private:
  
  CFuint m_nbDistanceZones;
  CFuint m_nbBasicFilters;
  CFreal m_inclusionCriteria;
}; // end of class InterpolationStencilComputer

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class defines a stencil of a basic filter used in interpolation-type
   * filters.
   *
   * @author Willem Deconinck
   */
class Triangle : public FilterStencil
{
private:
  
  /// Centroid of the triangle.
  CFuint m_centroid;
  
  /// Sub-area's of the triangle
  RealVector m_subArea;
  
  /// Total area
  CFreal m_totalArea;
  
  /// Coordinate Linker
  Common::SafePtr<CoordinateLinker> m_coordinateLinker; 

public:  
  
  /**
   * Default constructor
   * Should not be called! Necessary for std::map definitions
   */
  Triangle(): m_centroid(0), m_subArea(3) { /* should not be called */ }
  
  /**
   * Constructor
   */
  Triangle(const CFuint& centroid, Common::SafePtr<CoordinateLinker>& coordinateLinker ) :
    m_centroid(centroid), m_coordinateLinker(coordinateLinker), m_subArea(3) {}
    
  Triangle(const Triangle& T) : FilterStencil(T),
                                m_centroid(T.m_centroid),
                                m_subArea(T.m_subArea),
                                m_totalArea(T.m_totalArea),
                                m_coordinateLinker(T.m_coordinateLinker) {  }
    
  /**
   * Calculate the total area and the sub-area's
   * defined by the centroid en the 3 nodes of the triangle.
   */
  void calculateAreas() {
    CFreal s;
    CFreal lengthCentreTo1 = MathTools::MathFunctions::getDistance(m_coordinateLinker->getCoordinates(m_centroid),m_coordinateLinker->getCoordinates(this,0));
    CFreal lengthCentreTo2 = MathTools::MathFunctions::getDistance(m_coordinateLinker->getCoordinates(m_centroid),m_coordinateLinker->getCoordinates(this,0));
    CFreal lengthCentreTo3 = MathTools::MathFunctions::getDistance(m_coordinateLinker->getCoordinates(m_centroid),m_coordinateLinker->getCoordinates(this,0));
    CFreal length1To2 = MathTools::MathFunctions::getDistance(m_coordinateLinker->getCoordinates(this,0),m_coordinateLinker->getCoordinates(this,1));
    CFreal length2To3 = MathTools::MathFunctions::getDistance(m_coordinateLinker->getCoordinates(this,1),m_coordinateLinker->getCoordinates(this,2));
    CFreal length3To1 = MathTools::MathFunctions::getDistance(m_coordinateLinker->getCoordinates(this,2),m_coordinateLinker->getCoordinates(this,0));
    // largeTriangle
    s = 0.5*(length1To2 + length2To3 + length3To1);
    m_totalArea = sqrt(s*(s-length1To2)*(s-length2To3)*(s-length3To1));

    // subTriangle1
    s = 0.5*(lengthCentreTo1 + length1To2 + lengthCentreTo2);
    m_subArea[0] = sqrt(s*(s-lengthCentreTo1)*(s-length1To2)*(s-lengthCentreTo2));

    // subTriangle2
    s = 0.5*(lengthCentreTo2 + length2To3 + lengthCentreTo3);
    m_subArea[1] = sqrt(s*(s-lengthCentreTo2)*(s-length2To3)*(s-lengthCentreTo3));

    // subTriangle3
    s = 0.5*(lengthCentreTo3 + length3To1 + lengthCentreTo1);
    m_subArea[2] = sqrt(s*(s-lengthCentreTo3)*(s-length3To1)*(s-lengthCentreTo1));
  }

  /**
   * Calculate deviation of the ratio subArea/totalArea from 1/3
   * for all 3 subArea's and return the maximal deviation.
   * @return maximum deviation
   */
  CFreal calculateDeviation() {
    RealVector criteria(3);
    criteria = abs(3.*m_subArea/m_totalArea-1.);
    return criteria.max();
  }
  
  /**
   * @return The total area of the triangle
   */
  CFreal getTotalArea() {
    return m_totalArea;
  }
  
  /**
   * @return sub-area with index i
   */
  CFreal getSubArea(const CFuint& i){
    return m_subArea[i];
  }
  
  /**
   * @return The sum of the area's of the sub-triangles
   */
  CFreal calculateSumSubAreas() {
    return m_subArea.sum();
  }
  
  /**
   * @return The centroid to define the sub-triangles
   */
  CFuint getCentroid() {
    return m_centroid;
  }
  
  friend bool operator== (Triangle &T1, Triangle &T2);
  friend bool operator!= (Triangle &T1, Triangle &T2);
};

inline bool operator== (Triangle &T1, Triangle &T2) {
  return (T1.getElement(0) == T2.getElement(0) || T1.getElement(1) == T2.getElement(1) || T1.getElement(2) == T2.getElement(2));
}

inline bool operator!= (Triangle &T1, Triangle &T2) {
  return !(T1==T2);
}


//////////////////////////////////////////////////////////////////////////////

    } // end of namespace ExplicitFilters

  } // end of namespace Numerics

} // end of namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ExplicitFilters_InterpolationStencilComputer_hh
