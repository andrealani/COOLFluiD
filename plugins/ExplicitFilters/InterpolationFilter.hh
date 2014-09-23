#ifndef COOLFluiD_Numerics_ExplicitFilters_InterpolationFilter_hh
#define COOLFluiD_Numerics_ExplicitFilters_InterpolationFilter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/MethodStrategy.hh"
#include "Framework/State.hh"
#include "FilterStrategy.hh"
#include "InterpolationStencilComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
	
		namespace ExplicitFilters {
		  
//////////////////////////////////////////////////////////////////////////////

/** 
 * This class offers a basic interface for the different explicit
 * filter
 * 
 * @author Willem Deconinck
 */
class InterpolationFilter : public FilterStrategy {

public: // functions

	/** 
	 * Defines the Config Options of this class
	 * @param    options   an OptionList where to add the options
	 */
  static void defineConfigOptions(Config::OptionList& options);

  /** 
   * Returns the DataSocket's that this numerical strategy needs 
   * as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
    needsSockets();

  /** 
   * Constructor
   */
  InterpolationFilter(const std::string& name);

  /** 
   * Default destructor
   */
  virtual ~InterpolationFilter();

  /** 
   * Set private data that will be used during the computation
   */
  virtual void setup();

  /** 
   * Configure the object
   */
  virtual void configure(Config::ConfigArgs& args);

  /** 
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "InterpolationFilter";
  }
  
  virtual void getInfo();
  
protected: // helper functions
	
	
	/**
	 * Add a weight to the centre cell in the stencil, while
	 * preserving the filter properties, to ensure that the
	 * transfer function at the maximal wavenumber is zero.
	 * G(kmax) = 0
	 */
  void addCentreWeight(const CFuint& centreCellID);
  
  /**
   * Calculate the discrete filter weights
   * for a given cell
   * centreCellID The ID of the given cell
   */
  void calculateWeights();
  
public:

protected:
  
  CFreal calculateTriangleFilterMoment(const CFuint& q, const CFuint& r, Triangle& triangle, const RealVector& nodeWeights);
  
  
private: // data

  /// Flag for non negative least squares or regular least squares
  bool m_nnls;

	/// Factor to give weight to centre cell
	CFreal m_centreWeight;
	
	/// Factor to lift off transfer function
	CFreal m_relaxationFactor;
	
	/// Triangle stencils
	Framework::DataSocketSink < std::vector<Triangle> > socket_triangles;  

protected: // data


}; // end of class FilterStrategy

//////////////////////////////////////////////////////////////////////////////

  	} // end of namespace ExplicitFilters

	} // end of namespace Numerics
	
} // end of namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ExplicitFilters_InterpolationFilter_hh
