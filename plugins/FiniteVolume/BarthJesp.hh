#ifndef COOLFluiD_Numerics_FiniteVolume_BarthJesp_hh
#define COOLFluiD_Numerics_FiniteVolume_BarthJesp_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/Limiter.hh"
#include "Framework/Storage.hh"
#include "Framework/DataSocketSink.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

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
 * This class implements Barth and Jespersen limiter in 2D/3D for FVM
 *
 * @author Mehmet Sarp Yalim
 * @author Andrea Lani
 *
 */
class BarthJesp : public Framework::Limiter<CellCenterFVMData> {
public:
  
#ifdef CF_HAVE_CUDA
  /**
   * This nested class holds configurable options for this object
   *
   * @author Andrea Lani
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
    }
    
    CFreal alpha;
    bool useFullStencil;
  };
  
  /**
   * This nested class implements Barth and Jespersen limiter in 2D/3D for FVM
   * and can be also called inside GPU kernels
   *
   * @author Andrea Lani
   *
   */
  template <typename PHYS>
  class DeviceFunc {
  public:
    typedef BarthJesp BASE;
    
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
  } 
  
  /// copy the local configuration options to the Framework::DEVICE
  void copyConfigOptions(DeviceConfigOptions<NOTYPE>* dco) 
  {
    dco->useFullStencil = m_useFullStencil;
    dco->alpha = m_alpha;
  }
  
#endif
    
  /**
   * Constructor
   */
  BarthJesp(const std::string& name);
  
  /**
   * Default destructor
   */
  ~BarthJesp();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args )
  {
    Framework::Limiter<CellCenterFVMData>::configure(args);
  }

  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
  /**
   * Compute the limiter in the current face
   */
  void limit(const std::vector<std::vector<Framework::Node*> >& coord,
	     Framework::GeometricEntity* const cell,
	     CFreal* limiterValue);
  
private:
  
  /// storage for the stencil via pointers to neighbors
  Framework::DataSocketSink<std::vector<Framework::State*> > socket_stencil;
  
  /// socket for uX values
  Framework::DataSocketSink< CFreal> socket_uX;
  
  /// socket for uY values
  Framework::DataSocketSink< CFreal> socket_uY;
  
  /// socket for uZ values
  Framework::DataSocketSink< CFreal> socket_uZ;
  
}; // end of class BarthJesp

//////////////////////////////////////////////////////////////////////////////

#ifdef CF_HAVE_CUDA

template <typename PHYS>
void BarthJesp::DeviceFunc<PHYS>::limit(const KernelData<CFreal>* kd, 
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
  const CFreal* uZ = &kd->uZ[0];
  
  for(CFuint iVar = 0; iVar < PHYS::NBEQS; ++iVar) {
    CFreal min0 = state[iVar];
    CFreal max0 = state[iVar];
    
    // loop over the neighbor cells belonging to the chosen stencil
    const CFuint nbNeighbors = (m_dco->useFullStencil) ? stencilSize : nbFaces;
    for(CFuint in = 0; in < nbNeighbors; ++in) {
      const CFuint cellID = cell->getNeighborID(in);
      const CFint stype = cell->getNeighborType(in);
      const CFreal neighState = (stype > 0) ? states[cellID*PHYS::NBEQS+iVar] : 
    	ghostStates[cellID*PHYS::NBEQS+iVar];
      min0 = (min0 < neighState) ? min0 : neighState;
      max0 = (max0 > neighState) ? max0 : neighState;
    }
    
    const CFreal deltaPlusMax = max0 - state[iVar];
    const CFreal deltaPlusMin = min0 - state[iVar];
    CFreal psi = 1.0;
    CFreal psimin = 1.1;
    for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
      const CFuint startx = iFace*PHYS::DIM;
      const CFreal xq = coord[startx+XX];
      const CFreal yq = coord[startx+YY];
      const CFuint idx = stateID*PHYS::NBEQS + iVar;
      CFreal deltaMin = uX[idx]*(xq - stateCoord[XX]) + uY[idx]*(yq - stateCoord[YY]);
      if (PHYS::DIM == 3) {
    	const CFreal zq = coord[startx+ZZ];
    	deltaMin += uZ[idx]*(zq - stateCoord[ZZ]);
      }
      
      if (deltaMin > 0.0) {
    	const CFreal dMaxdMin = deltaPlusMax/deltaMin;
    	psi = (dMaxdMin < 1.) ? dMaxdMin : 1.0;
      }
      if (deltaMin < 0.0) {
    	const CFreal dMindMin = deltaPlusMin/deltaMin;
    	psi = (dMindMin < 1.) ? dMindMin : 1.0; 
      }
      psimin = (psi < psimin) ? psi : psimin;
    }
    
    limiterValue[iVar] = psimin;
  }
}
      
#endif   

//////////////////////////////////////////////////////////////////////////////
      
    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_BarthJesp_hh
