#ifndef COOLFluiD_Numerics_ExplicitFilters_FilterStrategy_hh
#define COOLFluiD_Numerics_ExplicitFilters_FilterStrategy_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/PhysicalModel.hh"
#include "Framework/BaseMethodStrategyProvider.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/MethodStrategy.hh"
#include "FilterData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
	
		namespace ExplicitFilters {
		        
      class FilterData;
      class FilterStencil;

//////////////////////////////////////////////////////////////////////////////

/** 
 * This class offers a basic interface for the different explicit
 * filter
 * 
 * @author Willem Deconinck
 */
class FilterStrategy : public Framework::MethodStrategy<FilterData> {
  friend class TransferFunctionWrapper;
  friend class Gfunc;
  friend class Gfunc_cutoff;

public: // functions

  enum direction {Xdirection,Ydirection,Zdirection,XYdiagonal,XZdiagonal,YZdiagonal,XYZdiagonal,Diagonal};

	/** 
	 * Defines the Config Options of this class
	 * @param    options   an OptionList where to add the options
	 */
  static void defineConfigOptions(Config::OptionList& options);

  typedef Framework::BaseMethodStrategyProvider
  <FilterData,FilterStrategy > PROVIDER;

  /** 
   * Constructor
   */
  FilterStrategy(const std::string& name);

  /** 
   * Default destructor
   */
  virtual ~FilterStrategy();

  /** 
   * Set private data that will be used during the computation
   */
  virtual void setup();

  /** 
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );

  /** 
   * Returns the DataSocket's that this numerical strategy needs 
   * as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
    needsSockets();
    
  /** 
   * Returns the DataSocket's that this numerical strategy needs 
   * as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> >
    providesSockets();

  /** 
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "FilterStrategy";
  }
  
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
  /**
   * Get info about the filter
   */
  virtual void getInfo()=0;

  /**
   * Output the transfer function
   */
  void outputTransferFunction();
  
  void calculateAllWeights();
  
  /**
   * Calculate the transfer function for a given stencil, and wave number
   * @param centreCellID The ID to identify the stencil
   * @param k The wave number
   * @return The complex value of the transfer function at k
   */
  CFcomplex transferFunction(RealVector& k);
  
  CFcomplex transferFunction(const RealVector& k, const std::vector<Framework::State*>& stencil, const RealVector& weights);

protected: // helper functions

  CFuint getNbWeightsToCompute();

  CFreal getFilterWidth();
  
  void setStencil(const Common::SafePtr<FilterStencil>& stencil)
  {
    m_stencil = stencil;
  }
  
  
  CFuint getStencilID() const
  {
    return m_stencil->getID();
  }
  
  Common::SafePtr<FilterStencil> getStencil() const
  {
    return m_stencil;
  }
  
  bool isGoodTransferfunction();

  CFreal targetTransferFunction(const RealVector& k);
  
  CFreal targetTransferFunction(const CFreal& kx, const CFreal& ky);
  
  CFreal targetTransferFunction(const CFreal& k);

  virtual void calculateWeights() = 0;
  
  virtual void outputAdditional() {}
  
  void calculateCentreCellWeight();

  RealVector getRealVectorInDirection(const direction dir);

  /**
   * Calculate maximal wavenumber visible on the given stencil
   * @param stencil The stencil or mesh
   * @return The maximal wavenumber
   */
  CFreal calculateMaximalWaveNumber();
  RealVector calculateMaximalWaveNumberVector();
  
  /**
   * Output the transfer function in gnuplot format
   * @param k_11 The wave numbers in xy-direction
   * @param G_11 The real transfer function values in xy-direction
   * @param k_10 The wave numbers in x-direction
   * @param G_10 The real transfer function values in x-direction
   * @param centreCellID The ID of the cell of the transfer function
   */
  void outputTransferFunctionGnuplot(RealVector& k_11, RealVector& G_11, RealVector& k_10, RealVector& G_10, CFuint& centreCellID);

  /**
   * Output the transfer function in tecplot format
   * @param K The wave numbers
   * @param G The complex transfer function values
   * @param N The number of samples in one direction
   * @param centreCellID The ID of the cell of the transfer function
   */
  void outputTransferFunctionTecplot(RealVector** K, CFcomplex** G, const CFuint& N, const CFuint& centreCellID) ;


  CFreal calculateFilterGridRatio(const direction dir);

  /**
   * Calculate the filter grid ratio of the given transfer function
   * based on the secon
   * Based on Lund (1997), On the use of discrete filters for large eddy simulation
   */
  CFreal calculateFilterGridRatioTransferFunctionMoment2();
  CFreal calculateFilterGridRatioTransferFunctionMoment20();
  CFreal calculateFilterGridRatioTransferFunctionMoment11();
  
  /**
   * Calculate the filter grid ratio of the given transfer function
   * based on the second moment
   * Based on Lund (1997), On the use of discrete filters for large eddy simulation
   */
  CFreal calculateFilterGridRatioFilterMoment2();
  
  /**
   * Calculate the filter moment
   */
  CFreal calculateFilterMoment(const CFuint& q, const CFuint& r);
  CFreal calculateFilterMoment(const CFuint& q, const CFuint& r, std::vector<Framework::State*>& stencil, RealVector& weights);
  
  CFreal calculateTransferFunctionMoment20();
  CFreal calculateTransferFunctionMoment20(const std::vector<Framework::State*>& stencil, const RealVector& weights);
  
  CFreal calculateTransferFunctionMoment11();
  CFreal calculateTransferFunctionMoment11(const std::vector<Framework::State*>& stencil, const RealVector& weights);
  
  
  CFreal calculateTransferFunctionMoment2();
  
  
  CFreal calculateFilterWaveNumberRootFinding();
  CFreal calculateFilterWaveNumberRootFinding(const direction dir);
  
  CFreal calculateMinimumTransferFunction();
  
  CFreal integrateTransferFunction(const RealVector &kmin, const RealVector& kmax);
  
  CFreal ky1(const CFreal& kx) { return 0.0; }
  
  CFreal ky2(const CFreal& kx) {
    CFreal kmax = calculateMaximalWaveNumber();
    return sqrt(kmax*kmax - kx*kx);
  }
  
private: // data

 	/// power of transfer function output
  CFreal m_transferFunctionPower;
  
  /// check if all weights have been calculated
  bool m_allWeightsCalculated;
  
  /// Allowed deviation from 0 of the transferfunction at the maximal wavenumber
  CFreal m_allowedDeviation;
  
  /// Allowed overshoot at half of the filter wavenumber
  CFreal m_allowedOvershoot;
  

protected: // data

  Common::SafePtr<FilterStencil> m_stencil;
  
  /// Temporary vector of weights, which will be copied into datastructure after its computations
  RealVector m_weights;

  /// temporary storage for the centre cell of the stencil
  Framework::State* m_centreCell;
    
  /// string that collects error messages during the weight calculation process
  std::string m_errorMessagesDuringWeightCalculation;
  
  /// Maximum number of recomputations of stencils and weights during weight calculation
  CFuint m_maxNbStencilRecomputations;
  
  /// Factor to lift off transfer function
	CFreal m_relaxationFactor;

}; // end of class FilterStrategy

//////////////////////////////////////////////////////////////////////////////

class Tol
{
public:
  bool operator() (const CFreal& fa, const CFreal& fb) {
    return (std::abs(fb-fa) < 0.001);
  }
};

//////////////////////////////////////////////////////////////////////////////

class Gfunc
{
  private:
    RealVector m_direction;
    FilterStrategy* m_filter;
  public:
    Gfunc(FilterStrategy* filter) : 
      m_filter(filter)
    { 
      CFuint dim = Framework::PhysicalModelStack::getActive()->getDim();
      m_direction.resize(dim);
      m_direction = m_filter->getRealVectorInDirection(FilterStrategy::XYdiagonal);
    }
  
    CFreal operator() (const CFreal &k) 
    {
      RealVector kvec = k*m_direction;
      return real(m_filter->transferFunction(kvec));
    }
};

//////////////////////////////////////////////////////////////////////////////

class Gfunc_cutoff
{
  private:
    RealVector m_direction;
    FilterStrategy* m_filter;
  public:
    Gfunc_cutoff(FilterStrategy* filter, const FilterStrategy::direction dir) : 
      m_filter(filter)
    { 
      CFuint dim = Framework::PhysicalModelStack::getActive()->getDim();
      m_direction.resize(dim);
      m_direction = m_filter->getRealVectorInDirection(dir);
    }
  
    CFreal operator() (const CFreal &k) 
    {
      RealVector kvec = k*m_direction;
      return real(m_filter->transferFunction(kvec)) - 0.5;
    }
};

//////////////////////////////////////////////////////////////////////////////

class TransferFunctionWrapper{
public:
    /* constructor */
    TransferFunctionWrapper(FilterStrategy *object,
                std::vector<Framework::State*>& stencil, 
                RealVector& weights)
      : Obj(object), m_stencil(stencil), m_weights(weights) { }
        
    // "member function evaluation" with one parameter
    CFreal operator() (const RealVector &k){
      return real((Obj->transferFunction)(k,m_stencil,m_weights));
    }
    CFreal operator() (const CFreal& kx, const CFreal& ky){
      RealVector k(2);
      k[XX] = kx;
      k[YY] = ky;
      return real((Obj->transferFunction)(k,m_stencil,m_weights));
    }
    
private:
    TransferFunctionWrapper();	/* make default constructor private */
    FilterStrategy *Obj;		/*!< pointer to the object */
    std::vector<Framework::State*>& m_stencil;
    RealVector& m_weights;
};


template<typename Object_Type, typename Member_Pointer>
class MemberFunctionWrapper{
public:
    /* constructor */
    MemberFunctionWrapper(Object_Type* object, Member_Pointer mem_func)
      : Obj(object), Ptr(mem_func) { }
        
    // "member function evaluation" with one parameter
    CFreal operator() (const CFreal& k){
      return (Obj->*Ptr)(k);
    }
    
private:
    // MemberFunctionWrapper(); /* make default constructor private */
    Object_Type *Obj;		/*!< pointer to the object */
    Member_Pointer Ptr;		  /*!< pointer to the class member function */
};

//////////////////////////////////////////////////////////////////////////////

  	} // end of namespace ExplicitFilters

	} // namespace Numerics
	
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ExplicitFilters_FilterStrategy_hh
