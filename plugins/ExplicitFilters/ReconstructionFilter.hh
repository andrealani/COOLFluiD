#ifndef COOLFluiD_Numerics_ExplicitFilters_ReconstructionFilter_hh
#define COOLFluiD_Numerics_ExplicitFilters_ReconstructionFilter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/MethodStrategy.hh"
#include "Framework/State.hh"
#include "FilterStrategy.hh"

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
class ReconstructionFilter : public FilterStrategy {
  friend class LHSWrapper;
  friend class RHSWrapper;
public: // functions

	/** 
	 * Defines the Config Options of this class
	 * @param    options   an OptionList where to add the options
	 */
  static void defineConfigOptions(Config::OptionList& options);

  /** 
   * Constructor
   */
  ReconstructionFilter(const std::string& name);

  /** 
   * Default destructor
   */
  virtual ~ReconstructionFilter();

  /** 
   * Set private data that will be used during the computation
   */
  virtual void setup();
  
  /** 
   * Returns the DataSocket's that this numerical strategy needs 
   * as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
    needsSockets();

  /** 
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );

  /** 
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "ReconstructionFilter";
  }
  
  virtual void getInfo();

  /**
   * Show the reconstruction equation
   *
   * @param   order   order of the reconstruction
   * @return  LaTeX formatted string for the equation used in the reconstruction
   */
  virtual std::string reconstructionEquation(const CFuint order) const =0;
  
	/**
	 * Calculate the number of unknowns in the least squares problem for
	 * a given order of polynomial reconstruction.
	 * @param order The order of reconstruction
	 * @return The number of unknowns
	 */
	static CFuint calculateNbUnknowns(const CFuint& order);

protected: // helper functions	
	
	/**
	 * Add a weight to the centre cell in the stencil, while
	 * preserving the filter properties, to ensure that the
	 * transfer function at the maximal wavenumber is zero.
	 * G(kmax) = 0
	 */
  void addRelaxationFactor();
  
  /**
   * Calculate the system matrix of the least squares problem Ax=b
   * @param A the system matrix to be calculated
   * @param stencil The stencil to be used for the calculation
   */
  virtual void calculateSystemMatrix(RealMatrix& A) =0;
  
  /**
   * Calculate the weights for weighted least squares reconstruction
   * @param W The diagonal matrix to store the weights in
   * @param stencil The stencil used to create the weights
   */
  void calculateWeightsMatrix(RealMatrix& W);
  
  /**
   * Calculate the discrete filter weights
   */
  void calculateAllWeights();
  
  /**
   * Calculate the discrete filter weights
   * for a given cell
   * centreCellID The ID of the given cell
   */
  void calculateWeights();
  
  void getFilterGridRatioVSweightingFactor();
  
  void leastSquaresReconstruction();

  void optimizedLeastSquaresReconstruction();
  
  virtual void outputAdditional();
  
private: // data
 
	/// Weighting factor for influencing filter width
	CFreal m_weightingFactor;
	
	/// Defines if weighting is used or not (default = true)
  bool m_weighting;
  
  /// Flag for non negative least squares or regular least squares
  bool m_nnls;

	/// Factor to give weight to centre cell
	CFreal m_centreWeight;
  
protected: // data


}; // end of class FilterStrategy

//////////////////////////////////////////////////////////////////////////////


class LHSWrapper{
public:
  /* constructor */
  LHSWrapper(ReconstructionFilter *object,
             const std::vector<Framework::State*>& stencil1, 
             const RealVector& weights1,
             const std::vector<Framework::State*>& stencil2,
             const RealVector& weights2)
    : Obj(object), m_stencil1(stencil1), m_weights1(weights1), m_stencil2(stencil2), m_weights2(weights2) { }
      
  // "member function evaluation" with one parameter
  CFreal operator() (const RealVector &k){
    return real((Obj->transferFunction)(k,m_stencil1,m_weights1))*real((Obj->transferFunction)(k,m_stencil2,m_weights2));
  }
  CFreal operator() (const CFreal &kx, const CFreal& ky){
    RealVector k(2);
    k[XX] = kx;
    k[YY] = ky;
    return real((Obj->transferFunction)(k,m_stencil1,m_weights1))*real((Obj->transferFunction)(k,m_stencil2,m_weights2));
  }
    
private:
  LHSWrapper();	/* make default constructor private */
  ReconstructionFilter *Obj;		/*!< pointer to the object */
  const std::vector<Framework::State*>& m_stencil1;
  const RealVector& m_weights1;
  const std::vector<Framework::State*>& m_stencil2;
  const RealVector& m_weights2;
};

class RHSWrapper{
public:
    /* constructor */
    RHSWrapper(ReconstructionFilter *object,
               const std::vector<Framework::State*>& stencil, 
               const RealVector& weights)
      : Obj(object), m_stencil(stencil), m_weights(weights) 
      {
        m_kmax = Obj->calculateMaximalWaveNumber();
      }
        
    // "member function evaluation" with one parameter
    CFreal operator() (const RealVector &k){
      return real((Obj->transferFunction)(k,m_stencil,m_weights)) * Obj->targetTransferFunction(k/m_kmax*MathTools::MathConsts::CFrealPi());
    }
    CFreal operator() (const CFreal &kx, const CFreal& ky){
      RealVector k(2);
      k[XX] = kx;
      k[YY] = ky;
      return real((Obj->transferFunction)(k,m_stencil,m_weights)) * Obj->targetTransferFunction(k/m_kmax*MathTools::MathConsts::CFrealPi());
    }
    
private:
    RHSWrapper();	/* make default constructor private */
    ReconstructionFilter *Obj;		/*!< pointer to the object */
    const std::vector<Framework::State*>& m_stencil;
    const RealVector& m_weights;
    CFreal m_kmax;
};

template<typename Member_Pointer>
class ReconstructionFilterMemberFunctionWrapper{
public:
    /* constructor */
    ReconstructionFilterMemberFunctionWrapper(ReconstructionFilter* object, Member_Pointer mem_func)
      : Obj(object), Ptr(mem_func) { }
        
    // "member function evaluation" with one parameter
    CFreal operator() (const CFreal& k){
      return (Obj->*Ptr)(k);
    }
    
private:
    // MemberFunctionWrapper(); /* make default constructor private */
    ReconstructionFilter *Obj;		/*!< pointer to the object */
    Member_Pointer Ptr;		  /*!< pointer to the class member function */
};

//////////////////////////////////////////////////////////////////////////////

  	} // namespace ExplicitFilters

	} // namespace Numerics
	
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_ExplicitFilters_ReconstructionFilter_hh
