#ifndef COOLFluiD_Numerics_FiniteVolume_Venktn2D_hh
#define COOLFluiD_Numerics_FiniteVolume_Venktn2D_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/Limiter.hh"
#include "Framework/Storage.hh"
#include "Framework/DataSocketSink.hh"
#include "CellCenterFVMData.hh"

#ifdef CF_HAVE_CUDA
#include "FiniteVolume/CellData.hh"
#include "FiniteVolume/FluxData.hh"
#include "FiniteVolume/KernelData.hh"
#include "Common/CUDA/CudaEnv.hh"
#include "Common/CUDA/CFVec.hh"
#endif

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class GeometricEntity;
  }

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class implements Venkatakrishnan limiter in 2D for FVM
 *
 * @author Mehmet Sarp Yalim
 * @author Isaac Alonso Asensio (GPU-port)
 */
class Venktn2D : public Framework::Limiter<CellCenterFVMData> {
public:

#ifdef CF_HAVE_CUDA
  /**
   * This nested class holds configurable options for this object
   *
   * @author Isaac Alonso
   *
   */
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
      alpha = in->alpha;
      useFullStencil = in->useFullStencil;
      deltaMin = in->deltaMin;
      for (CFuint i=0; i<18; i++) {
         magnitudeValues[i] = in->magnitudeValues[i];
      }
      coeffEps = in->coeffEps;
      length = in->length;
    }
    
    CFreal alpha;
    bool useFullStencil;
    CFreal deltaMin;
    CFreal magnitudeValues[18]; //IA: Change for P::NBEQS and when defining dco @ FVM-CUDA include VarSet as template parameter
    CFreal coeffEps;
    CFreal length;
    bool isMFMHD; 
  };
  
  /**
   * This nested class implements Venkatakrishnan limiter in 2D for FVM
   * and can be also called inside GPU kernels
   *
   * @author Isaac Alonso
   *
   */
  template <typename PHYS>
  class DeviceFunc {
  public:
    typedef Venktn2D BASE;
    
    /// constructor taking the options as argument
    HOST_DEVICE DeviceFunc(DeviceConfigOptions<NOTYPE>* dco) : m_dco(dco) {}
    
    /// compute the flux limiter values (one per variables)
    HOST_DEVICE void limit(const KernelData<CFreal>* kd, 
			   const CellData::Itr* cell,
			   const CFreal* coord, 
			   CFreal* limiterValue);
    
  private:
    /// options to be used on the Framework::DEVICE
    DeviceConfigOptions<NOTYPE>* m_dco;
  };
  
  /// copy the local configuration options to the Framework::DEVICE
  void copyConfigOptionsToDevice(DeviceConfigOptions<NOTYPE>* dco) 
  {
   
    CudaEnv::copyHost2Dev(&dco->alpha, &m_alpha, 1);
    CudaEnv::copyHost2Dev(&dco->useFullStencil, &m_useFullStencil, 1);
    CudaEnv::copyHost2Dev(&dco->deltaMin, &_deltaMin, 1);
    CudaEnv::copyHost2Dev(&dco->magnitudeValues[0], &_magnitudeValues[0], 18);
    CudaEnv::copyHost2Dev(&dco->coeffEps, &_coeffEps, 1);
    CudaEnv::copyHost2Dev(&dco->length, &_length, 1);
    CudaEnv::copyHost2Dev(&dco->isMFMHD, &_isMFMHD, 1);
  } 
  
  /// copy the local configuration options to the Framework::DEVICE
  void copyConfigOptions(DeviceConfigOptions<NOTYPE>* dco) 
  {
    dco->alpha = m_alpha;
    dco->useFullStencil = m_useFullStencil;
    dco->deltaMin = _deltaMin;
    for (CFuint i=0; i<18; i++) {
       dco->magnitudeValues[i] = _magnitudeValues[i];
    }
    dco->coeffEps = _coeffEps;
    dco->length = _length;
    dco->isMFMHD = _isMFMHD;
  }
  
#endif


  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  Venktn2D(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~Venktn2D();
  
  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets()
  {
    std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = 
      Framework::Limiter<CellCenterFVMData>::needsSockets();

    result.push_back(&socket_stencil);
    result.push_back(&socket_uX);
    result.push_back(&socket_uY);
    return result;
  }

  /**
   * Compute the limiter in the current face
   */
  virtual void limit(const std::vector<std::vector<Framework::Node*> >& coord,
		     Framework::GeometricEntity* const cell,
		     CFreal* limiterValue);

  /**
   * Set up the private data
   */
  virtual void setup();
  
protected:
  
  /**
   * Compute the denominator of the limiter argument
   */
  void computeDeltaMin(const Framework::Node& coord, const Framework::State& state, CFuint iVar)
  {
    const CFuint stateID = state.getLocalID();
    const RealVector& stateCoord = state.getCoordinates();
    _deltaMin = (socket_uX.getDataHandle())(stateID,iVar,state.size())*(coord[XX] - stateCoord[XX]) + 
      (socket_uY.getDataHandle())(stateID,iVar,state.size())*(coord[YY] - stateCoord[YY]);
  }
  
protected:

  /// storage for the stencil via pointers to neighbors
  Framework::DataSocketSink<std::vector<Framework::State*> > socket_stencil;
  
  /// socket for uX values
  Framework::DataSocketSink<CFreal> socket_uX;
  
  /// socket for uY values
  Framework::DataSocketSink<CFreal> socket_uY;
  
  /// denominator of limiter argument 
  CFreal _deltaMin;
  
  /// user defined state vector with order of magnitude of the solution 
  std::vector<CFreal> _magnitudeValues;
  
  /// user defined coefficient for the epsilon
  CFreal _coeffEps;

  /// user defined characteristic solution length in the smooth flow region
  CFreal _length;
 
  /// Flag for MultiFluid cases to use the fix for smooth regions 
  bool _isMFMHD; 
}; // end of class Venktn2D

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_CUDA

template <typename PHYS>
void Venktn2D::DeviceFunc<PHYS>::limit(const KernelData<CFreal>* kd, 
				       const CellData::Itr* cell,
				       const CFreal* coord,
				       CFreal* limiterValue)
{
  using namespace COOLFluiD::CudaEnv;
  
  const CFuint stencilSize = cell->getStencilSize();
  const CFuint nbFaces = cell->getNbActiveFacesInCell();
  const CFuint stateID = cell->getCellID();
  const CFreal* state = &kd->states[stateID*PHYS::NBEQS];
  const CFreal* stateCoord = &kd->centerNodes[stateID*PHYS::DIM];
  const CFreal* states = &kd->states[0];
  const CFreal* ghostStates = &kd->ghostStates[0];
  const CFreal* uX = &kd->uX[0];
  const CFreal* uY = &kd->uY[0];
  //const CFreal* uZ = &kd->uZ[0];
   
//  if (stateID == 319) { printf("StateCoord \t XX %.15e \t YY %.15e \n", stateCoord[0], stateCoord[1]); }
  for(CFuint iVar = 0; iVar < PHYS::NBEQS; ++iVar) {
    CFreal min0 = state[iVar];
    CFreal max0 = state[iVar];
    CFreal avgDistance = 0.0;
    
    // loop over the neighbor cells belonging to the chosen stencil
    const CFuint nbNeighbors = (m_dco->useFullStencil) ? stencilSize : nbFaces;
    for(CFuint in = 0; in < nbNeighbors; ++in) {
      const CFuint cellID = cell->getNeighborID(in);
      const CFint stype = cell->getNeighborType(in);
      const CFreal neighState = (stype > 0) ? states[cellID*PHYS::NBEQS+iVar] : 
    	ghostStates[cellID*PHYS::NBEQS+iVar];
      const CFreal* neighCoord = (stype > 0) ? &kd->centerNodes[cellID*PHYS::DIM] :
        &kd->ghostNodes[cellID*PHYS::DIM];

      min0 = (min0 < neighState) ? min0 : neighState;
      max0 = (max0 > neighState) ? max0 : neighState;
      //avgDistance += MathFunctions::getDistance(state->getCoordinates(), currState->getCoordinates());
      // getDistance() is not a CUDA adapted function. So it needs to be reimplemented
      CFreal dist = 0.;
      for (CFuint i = 0; i < PHYS::DIM; ++i)
      {
         const CFreal diff = neighCoord[i] - stateCoord[i];
         dist += diff*diff;
      }
//      if (stateID == 319) { printf("NeighCoord \t XX %.15e \t YY %.15e \n", neighCoord[0], neighCoord[1]); }
//      if (stateID == 319) { printf("NeighCoord \t XX %.15e \t YY %.15e @ Cell %d\n", neighCoord[0], neighCoord[1], cellID); }
      avgDistance += sqrt(dist);

    }
    avgDistance /= nbNeighbors;
   
    const CFreal deltaPlusMax = max0 - state[iVar];
    const CFreal deltaPlusMin = min0 - state[iVar];
   
    
 
    CFreal psi = 1.0;
    CFreal psimin = 1.1;
    if(m_dco->isMFMHD){ // IMPLEMENTATION FROM VENKATAKRISHNAN PAPER 
      const CFreal epsilonFactor = m_dco->coeffEps*avgDistance/m_dco->length;
      const CFreal epsilon2 = epsilonFactor*epsilonFactor*epsilonFactor; //_magnitudeValues[iVar]*_magnitudeValues[iVar]*
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
        // number of quadrature points associated with this face
        const CFuint nbPoints = 1;
      
        for (CFuint ip = 0; ip < nbPoints; ++ip) {
          const CFuint startx = iFace*PHYS::DIM;
          const CFreal xq = coord[startx+XX];
          const CFreal yq = coord[startx+YY];
          const CFuint idx = stateID*PHYS::NBEQS + iVar;
          CFreal deltaMin = uX[idx]*(xq - stateCoord[XX]) + uY[idx]*(yq - stateCoord[YY]);
         
          /*  IA: TODO, check this implementation, which could be faster as it avoid 'if' statements!!
          const bool sign = (deltaMin < 0) ?  false : true;   // True if deltamin is positive
          const CFreal dPlus = (sign) ? deltaPlusMax/m_dco->magnitudeValues[iVar] : deltaPlusMin/m_dco->magnitudeValues[iVar];	
          const CFreal dPlus2 = dPlus*dPlus;
          const CFreal dMinStar = (sign) ? (deltaMin/m_dco->magnitudeValues[iVar] + 1e-24) : -(deltaMin/m_dco->magnitudeValues[iVar]);
	  const CFreal dPlusMin = dPlus*dMinStar;
          const CFreal Num = (dPlus2 + epsilon2)*dMinStar + 2*dMinStar*dMinStar*dPlus;
          const CFreal Den = (dPlus2 + 2*dMinStar*dMinStar + dPlusMin + epsilon2);
          psi = 1./dMinStar*(Num/Den);
          printf("psi[%d] = %.15e \n", iVar, psi);
          */
          
          if (deltaMin > 0.0) {
	    const CFreal dMinStar = (fabs(deltaMin)/m_dco->magnitudeValues[iVar] + 1e-24);
	    const CFreal dPlus = deltaPlusMax/m_dco->magnitudeValues[iVar];
	    const CFreal dPlus2   = dPlus*dPlus;
	    const CFreal dPlusMin = dPlus*dMinStar;
	    const CFreal Num = (dPlus2 + epsilon2)*dMinStar + 2*dMinStar*dMinStar*dPlus;
	    const CFreal Den = (dPlus2 + 2*dMinStar*dMinStar + dPlusMin + epsilon2);
	    psi = 1./dMinStar*(Num/Den);
	  }
          if (deltaMin < 0.0) {
	    const CFreal dMinStar = -(fabs(deltaMin)/m_dco->magnitudeValues[iVar] + 1e-24);
            const CFreal dPlus = deltaPlusMin/m_dco->magnitudeValues[iVar];
            const CFreal dPlus2   = dPlus*dPlus;
            const CFreal dPlusMin = dPlus*dMinStar;
            const CFreal Num = (dPlus2 + epsilon2)*dMinStar + 2*dMinStar*dMinStar*dPlus;
            const CFreal Den = (dPlus2 + 2*dMinStar*dMinStar + dPlusMin + epsilon2);
            psi = 1./dMinStar*(Num/Den);
	  }
          
          /* // DEBUG
          printf("nbNeighbors = %d ; stencilSize = %d ; nbFaces = %d \n", nbNeighbors, stencilSize, nbFaces);
          printf("epsilon2 = %.15e ; _coeffEps = %.15e ; avgDistance = %.15e _length = %.15e \n", epsilon2, m_dco->coeffEps, avgDistance, m_dco->length); 
          printf("deltaPlusMax  %.15e deltaPlusMin %.15e \n", deltaPlusMax, deltaPlusMin); 
          printf("deltaMin[%d] = %.15e \n", iVar, deltaMin); 
          printf("psi[%d] = %.15e \n", iVar, psi); 
          */
          psimin = (psimin<psi) ? psimin :  psi;
        }
      }
    }
    /* IA: TODO, for no multifluidMHD cases
    else{ //OLD IMPLEMENTATION
      const CFreal epsilon = _coeffEps*_magnitudeValues[iVar]*_magnitudeValues[iVar]*
        pow(avgDistance/_length, 3.0);

      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
        const CFuint nbPoints = 1;

        for (CFuint ip = 0; ip < nbPoints; ++ip) {
          computeDeltaMin((*coord[iFace][ip]), *state, iVar);
 
          if (_deltaMin > 0.0) {
            const CFreal dPlusMax2   = deltaPlusMax*deltaPlusMax;
            const CFreal dPlusMaxMin = deltaPlusMax*_deltaMin;
            psi = (dPlusMax2 + 2.0*dPlusMaxMin + epsilon)/
              (dPlusMax2 + dPlusMaxMin + 2.0*_deltaMin*_deltaMin + epsilon);
          }
          if (_deltaMin < 0.0) {
            const CFreal y = deltaPlusMin/_deltaMin;
            psi = (y*y + 2*y)/(y*y + y + 2.);
          }
          psimin = min(psi, psimin);
        }
      }
    }
    */
    limiterValue[iVar] = psimin;

    
    const CFreal maxAllowableLimiterFunctionValue = 1.094;
    /* // DEBUG
    if (stateID == 319) { printf("limiterValue[%d] = %.15e @ Cell %d \n", iVar, psimin, stateID); }
    */
    if (limiterValue[iVar] > maxAllowableLimiterFunctionValue) {
      printf("wrong limiterValue %.15e \n", limiterValue[iVar]);
    }
  }
}
      
#endif   

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_Venktn2D_hh
