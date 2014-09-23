#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/SubOutletEuler2DPuvtPis.hh"

#include <iostream>
#include <cmath>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SubOutletEuler2DPuvtPis, CellCenterFVMData, FiniteVolumeNavierStokesModule>
subOutletEuler2DPuvtPisFVMCCProvider("SubOutletEuler2DPuvtPisFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler2DPuvtPis::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal > ("Pis","P_is");
}

//////////////////////////////////////////////////////////////////////////////

SubOutletEuler2DPuvtPis::SubOutletEuler2DPuvtPis(const string& name) :
  FVMCC_BC(name)
{
  addConfigOptionsTo(this);
  P_is = 0.0;
  setParameter("Pis",&P_is);
}

//////////////////////////////////////////////////////////////////////////////

SubOutletEuler2DPuvtPis::~SubOutletEuler2DPuvtPis()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler2DPuvtPis::setGhostState(GeometricEntity *const face)
{

  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  
  const CFuint faceID = face->getID();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  CFreal nx = normals[startID];
  CFreal ny = normals[startID + 1];
  const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny);
  nx *= invFaceLength;
  ny *= invFaceLength;
  const CFreal u = (*innerState)[1];
  const CFreal v = (*innerState)[2];
  const CFreal gamma = (getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>())->getModel()->getGamma();
  const CFreal R = (getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>())->getModel()->getR();
  const CFreal aInnerState = sqrt(gamma*R*(*innerState)[3]);
  const CFreal vn = u*nx + v*ny;
  const CFreal Vel = sqrt(u*u + v*v);
  const CFreal M_normal_isoentropic = vn / aInnerState;
         
  //subsonic outlet face case 
  if(M_normal_isoentropic<1.0){
    //find the average pressure on the cells along the outlet boundary
    FVMCC_BC::setup();
    SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
    Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> _faceBuilder;
    _faceBuilder.setup();
    _faceBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
    FaceTrsGeoBuilder::GeoData& faceData = _faceBuilder.getDataGE();
    faceData.isBFace = true;
    faceData.trs = trs;
    const CFuint nbTrsFaces = trs->getLocalNbGeoEnts();
    CFreal P_sum = 0.0;  
    CFreal T_sum = 0.0;
    CFreal M_sum = 0.0;
    
    for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace){
      CFLogDebugMed( "iFace = " << iFace << "\n");
      // build the GeometricEntity
      faceData.idx = iFace;
      GeometricEntity *const iface = _faceBuilder.buildGE();
      const CFuint faceGlobalID = iface->getID();
      State *const InnerState = iface->getState(0);
      P_sum = P_sum + (*InnerState)[0];
      T_sum = T_sum + (*InnerState)[3];
      M_sum = M_sum + sqrt(((*InnerState)[1]*(*InnerState)[1]+(*InnerState)[2]*(*InnerState)[2])/(gamma*R*(*InnerState)[3]));
      // release the GeometricEntity
      _faceBuilder.releaseGE();
    }
    const CFreal P_average_inner = P_sum/nbTrsFaces;
    const CFreal T_average_inner = T_sum/nbTrsFaces;
    const CFreal M_average_inner = M_sum/nbTrsFaces;
    CFreal ratio_ave = P_average_inner/P_is;
    if(ratio_ave>1){
      ratio_ave = 1 / ratio_ave;
    }
    const CFreal damping_factor = 1;
    const CFreal damping = pow((1 - abs(ratio_ave - 1)),2);
    const CFreal P_exit = P_is + ((*innerState)[0] - P_average_inner)*damping;
    const CFreal P_ghost = 2*P_exit - (*innerState)[0];
    std::cout <<" Convergence condition on the outlet:  " <<std::endl;
    std::cout <<" P_is:  " <<P_is  <<" P_average_inner:  " <<P_average_inner  <<" ratio_ave:  " <<ratio_ave  <<" damping:  " <<damping  <<std::endl;

    // set the variables in the ghost cell
    (*ghostState)[0] = P_ghost;
    (*ghostState)[1] = (*innerState)[1];
    (*ghostState)[2] = (*innerState)[2];
    (*ghostState)[3] = (*innerState)[3];
  }
  //supersonic outlet face case
  if(M_normal_isoentropic>=1.0){
    (*ghostState)[0] = (*innerState)[0];
    (*ghostState)[1] = (*innerState)[1];
    (*ghostState)[2] = (*innerState)[2];
    (*ghostState)[3] = (*innerState)[3];
  }

}

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler2DPuvtPis::setup()
{

  FVMCC_BC::setup();
  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
