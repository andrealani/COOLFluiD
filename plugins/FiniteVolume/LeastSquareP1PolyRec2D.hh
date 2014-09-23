#ifndef COOLFluiD_Numerics_FiniteVolume_LeastSquareP1PolyRec2D_hh
#define COOLFluiD_Numerics_FiniteVolume_LeastSquareP1PolyRec2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolume/FVMCC_PolyRec.hh"

#ifdef CF_HAVE_CUDA
#include "FiniteVolume/FluxData.hh"
#include "FiniteVolume/KernelData.hh"
#include "FiniteVolume/CellData.hh"
#include "Common/CUDA/CudaEnv.hh"
#include "Framework/SubSystemStatus.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class implements a least square polynomial reconstructor in 2D for FVM
 *
 * @author Mehmet Sarp Yalim
 * @author Andrea Lani
 */
class LeastSquareP1PolyRec2D : public FVMCC_PolyRec {
public:

#ifdef CF_HAVE_CUDA
  /// nested class defining local options
  template <typename P = NOTYPE>
  class DeviceConfigOptions {
  public:
    /// constructor
    HOST_DEVICE DeviceConfigOptions() {}
    
    /// destructor
    HOST_DEVICE ~DeviceConfigOptions() {}
    
    /// initialize with another object of the same kind
    HOST_DEVICE void init(DeviceConfigOptions<P> *const in) 
    {
      limitIter      = in->limitIter;
      currIter       = in->currIter; 
      limitRes       = in->limitRes; 
      currRes        = in->currRes;  
      gradientFactor = in->gradientFactor;
      freezeLimiter  = in->freezeLimiter;
    }
    
    CFuint limitIter;      /// iteration at which limiter is activated
    CFuint currIter;       /// current iteration
    CFreal limitRes;       /// residual at which limiter is activated
    CFreal currRes;        /// current residual
    CFreal gradientFactor; /// gradient factor within [0,1]
    CFuint freezeLimiter;    /// flag controlling the freezing of the limiter
  };
  
  /// nested class defining a functor
  template <typename PHYS>
  class DeviceFunc {
  public:
    typedef LeastSquareP1PolyRec2D BASE;

    /// constructor taking the options as argument
    HOST_DEVICE DeviceFunc(DeviceConfigOptions<NOTYPE>* dco) : m_dco(dco) {}
    
    /// compute the cell-based gradients
    HOST_DEVICE void computeGradients(const CFreal* state, const CFreal* node,
				      KernelData<CFreal> *const kd,
				      CellData::Itr* cell); 
    
    /// extrapolate the solution on the given face
    HOST_DEVICE void extrapolateOnFace(FluxData<PHYS>* currFd, const CFreal* coord, const CFreal* uX, 
				       const CFreal* uY, const CFreal* uZ, const CFreal* limiter);
    
    /// extrapolate the solution on the given face for a given variable
    HOST_DEVICE inline void extrapolateOnFace(const CFuint iVar, 
					      FluxData<PHYS>* currFd, const CFreal* coord, const CFreal* uX, 
					      const CFreal* uY, const CFreal* uZ, const CFreal* limiter);
    
    /// compute the distance between two points
    template <CFuint DIM, typename T1, typename T2>
    HOST_DEVICE inline CFreal getDistance(const T1* n1, const T2* n2)
    {
      CFreal dist = 0.;
      for (size_t i = 0; i < DIM; ++i) {
	const CFreal diff = n1[i] - n2[i]; dist += diff*diff;
      }
      return sqrt(dist);
    }
    
  private:
    DeviceConfigOptions<NOTYPE>* m_dco;
  };
  
  /// copy the local configuration options to the device
  void copyConfigOptionsToDevice(DeviceConfigOptions<NOTYPE>* dco) 
  {
    // consider to copy to constant memory
    CudaEnv::copyHost2Dev(&dco->limitIter, &_limitIter, 1);
    CudaEnv::copyHost2Dev(&dco->limitRes, &_limitRes, 1);
    CudaEnv::copyHost2Dev(&dco->gradientFactor, &_gradientFactor, 1);
    
    CFuint flag = (_freezeLimiter) ? 1 : 0;
    CudaEnv::copyHost2Dev(&dco->freezeLimiter, &flag, 1);
    // global parameters for monitoring the simulation
    CFuint iter = Framework::SubSystemStatusStack::getActive()->getNbIter();
    CFreal res  = Framework::SubSystemStatusStack::getActive()->getResidual();    
    CudaEnv::copyHost2Dev(&dco->currIter, &iter, 1);
    CudaEnv::copyHost2Dev(&dco->currRes, &res, 1);
  }   
  
  /// copy the local configuration options to the device
  void copyConfigOptions(DeviceConfigOptions<NOTYPE>* dco) 
  {
    // consider to copy to constant memory
    dco->limitIter = _limitIter;
    dco->limitRes = _limitRes;
    dco->gradientFactor = _gradientFactor;
    dco->freezeLimiter = static_cast<CFuint>(_freezeLimiter);
    // global parameters for monitoring the simulation
    dco->currIter = Framework::SubSystemStatusStack::getActive()->getNbIter();
    dco->currRes = Framework::SubSystemStatusStack::getActive()->getResidual();
  }   
#endif
  
  /**
   * Constructor
   */
  LeastSquareP1PolyRec2D(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~LeastSquareP1PolyRec2D();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
  needsSockets();

  /**
   * Compute the gradients
   */
  virtual void computeGradients();

  /**
   * Set up the private data
   */
  virtual void setup();

  /**
   * Update the weights when nodes are moving
   */
  virtual void updateWeights();

protected:

  /**
   * Extrapolate the solution in the face quadrature points
   */
  virtual void extrapolateImpl(Framework::GeometricEntity* const face);

  /**
   * Extrapolate the solution in the face quadrature points
   */
  virtual void extrapolateImpl(Framework::GeometricEntity* const face,
			       CFuint iVar, CFuint leftOrRight);
  
protected:

  /// socket for stencil
  Framework::DataSocketSink<std::vector<Framework::State*> > socket_stencil;

  /// socket for weights
  Framework::DataSocketSink<CFreal> socket_weights;

  /// socket for uX values
  Framework::DataSocketSink<CFreal> socket_uX;

  /// socket for uY values
  Framework::DataSocketSink<CFreal> socket_uY;

  RealVector  _l11;

  RealVector  _l12;

  RealVector  _l22;

  RealVector  _lf1;

  RealVector  _lf2;

}; // end of class LeastSquareP1PolyRec2D

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_CUDA

template <typename PHYS>
void LeastSquareP1PolyRec2D::DeviceFunc<PHYS>::computeGradients
(const CFreal* state,  const CFreal* node, KernelData<CFreal> *const kd, CellData::Itr* cell)
{
  CFreal l11 = 0.0;
  CFreal l12 = 0.0;
  CFreal l22 = 0.0;
  CFreal lf1[PHYS::NBEQS];
  CFreal lf2[PHYS::NBEQS];
  for (CFuint i = 0; i< PHYS::NBEQS; ++i) {
    lf1[i] = 0.;
    lf2[i] = 0.;
  }  
  
  const CFreal* states = &kd->states[0];
  const CFreal* ghostStates = &kd->ghostStates[0];
  const CFreal* centerNodes = &kd->centerNodes[0];
  const CFreal* ghostNodes = &kd->ghostNodes[0];
  
  // loop over the neighbor cells belonging to the chosen stencil
  const CFuint stencilSize = cell->getStencilSize();
  for(CFuint in = 0; in < stencilSize; ++in) {
    const CFuint cellID = cell->getNeighborID(in);
    const CFint stype = cell->getNeighborType(in);
    // here a check on stype != 0 should not be needed since all
    // partition faces are neglected when stencil is constructed
    const CFreal* nodeLast = (stype > 0) ? &centerNodes[cellID*PHYS::DIM] : &ghostNodes[cellID*PHYS::DIM];  
    const CFreal deltaR = getDistance<PHYS::DIM, CFreal, CFreal>(node, nodeLast);
    const CFreal weight = 1.0/deltaR;
    const CFreal dx = weight*(nodeLast[XX] - node[XX]);
    const CFreal dy = weight*(nodeLast[YY] - node[YY]);
    l11 += dx*dx;
    l12 += dx*dy;
    l22 += dy*dy;
    
    const CFreal* neighState = (stype > 0) ? &states[cellID*PHYS::NBEQS] : &ghostStates[cellID*PHYS::NBEQS];  
    for (CFuint iVar = 0; iVar < PHYS::NBEQS; ++iVar) {
      const CFreal du = weight*(neighState[iVar] - state[iVar]);
      lf1[iVar] += dx*du;
      lf2[iVar] += dy*du;
    }
  }
  
  const CFreal det = l11*l22 - l12*l12;
  const CFreal invDet = 1./det; 
  CFreal* uX = &kd->uX[0];
  CFreal* uY = &kd->uY[0];
  const CFuint starts = cell->getCellID()*PHYS::NBEQS;
  for (CFuint i = 0; i < PHYS::NBEQS; ++i) {
    const CFuint gradx = starts + i;
    uX[gradx] = (l22*lf1[i] - l12*lf2[i])*invDet;
    uY[gradx] = (l11*lf2[i] - l12*lf1[i])*invDet;
  }
}
      
//////////////////////////////////////////////////////////////////////////////

template <typename PHYS>
void LeastSquareP1PolyRec2D::DeviceFunc<PHYS>::extrapolateOnFace
(FluxData<PHYS>* currFd, const CFreal* coord, const CFreal* uX, 
 const CFreal* uY, const CFreal* uZ, const CFreal* limiter)
{
  // gradient is not computed for ghost states (the corresponding inner gradient is taken instead)
  const CFuint stateIDL = currFd->getStateID(LEFT);
  const CFuint stateIDR = (!currFd->isBFace()) ? currFd->getStateID(RIGHT) : stateIDL;
  const CFreal* stateR = currFd->getState(RIGHT);
  const CFreal* stateL = currFd->getState(LEFT);
  CFreal* rstateR = currFd->getRstate(RIGHT);
  CFreal* rstateL = currFd->getRstate(LEFT);
  CFreal* nodeR = currFd->getNode(RIGHT);
  CFreal* nodeL = currFd->getNode(LEFT);
  
  // L/R cell coordinates are overwritten with the quadrature point 
  CFreal* rnodeR = currFd->getRnode(RIGHT);
  CFreal* rnodeL = currFd->getRnode(LEFT);
  for (CFuint i = 0; i < PHYS::DIM; ++i) {
    rnodeR[i] = rnodeL[i] = coord[i];
  }
  
  // AL: the following is more consistent in first order
  // for (CFuint i = 0; i < PHYS::DIM; ++i) {
  //   rnodeR[i] = (m_dco->gradientFactor > 0.) ? coord[i] : nodeR[i];
  //   rnodeL[i] = (m_dco->gradientFactor > 0.) ? coord[i] : nodeL[i];
  // }
  
  const CFuint startR = stateIDR*PHYS::NBEQS;
  const CFuint startL = stateIDL*PHYS::NBEQS;
  for (CFuint iVar = 0; iVar < PHYS::NBEQS; ++iVar) {
    const CFuint sR = startR + iVar;
    const CFreal duxR = uX[sR]*(coord[XX] - nodeR[XX]);
    const CFreal duyR = uY[sR]*(coord[YY] - nodeR[YY]);
    const CFreal duTotR = duxR + duyR;
    rstateR[iVar] = stateR[iVar];
    rstateR[iVar] += m_dco->gradientFactor*limiter[sR]*duTotR;
    const CFuint sL = startL + iVar;
    const CFreal duxL = uX[sL]*(coord[XX] - nodeL[XX]);
    const CFreal duyL = uY[sL]*(coord[YY] - nodeL[YY]);
    const CFreal duTotL = duxL + duyL;
    rstateL[iVar] = stateL[iVar];
    rstateL[iVar] += m_dco->gradientFactor*limiter[sL]*duTotL;
  }
}

//////////////////////////////////////////////////////////////////////////////

template <typename PHYS>
inline void LeastSquareP1PolyRec2D::DeviceFunc<PHYS>::extrapolateOnFace
(const CFuint iVar, FluxData<PHYS>* currFd, const CFreal* coord, const CFreal* uX, 
 const CFreal* uY, const CFreal* uZ, const CFreal* limiter)
{
  // gradient is not computed for ghost states (the corresponding inner gradient is taken instead)
  const CFuint sL   = currFd->getStateID(LEFT)*PHYS::NBEQS + iVar;
  const CFreal* nodeL = currFd->getNode(LEFT);
  const CFreal duxL = uX[sL]*(coord[XX] - nodeL[XX]);
  const CFreal duyL = uY[sL]*(coord[YY] - nodeL[YY]);
  const CFreal duTotL = duxL + duyL;
  currFd->getRstate(LEFT)[iVar] = currFd->getState(LEFT)[iVar];
  currFd->getRstate(LEFT)[iVar] += m_dco->gradientFactor*limiter[sL]*duTotL;
}

#endif
      
//////////////////////////////////////////////////////////////////////////////
      
    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_LeastSquareP1PolyRec2D_hh
