#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/SubOutletEuler2DPuvtMisMPI.hh"

#include <iostream>
#include <cmath>

#include "mpi.h"
#include "Common/PE.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SubOutletEuler2DPuvtMisMPI, CellCenterFVMData, FiniteVolumeNavierStokesModule>
subOutletEuler2DPuvtMisMPIFVMCCProvider("SubOutletEuler2DPuvtMisMPIFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler2DPuvtMisMPI::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal > ("Mis","M_is");
  options.addConfigOption< CFreal > ("Pt_inlet","Pt_inlet");
}

//////////////////////////////////////////////////////////////////////////////

SubOutletEuler2DPuvtMisMPI::SubOutletEuler2DPuvtMisMPI(const std::string& name) :
  FVMCC_BC(name),
  _comm(),
  _myP(),
  _nP()
{
  addConfigOptionsTo(this);
  M_is = 0.0;
  setParameter("Mis",&M_is);
  Pt_inlet = 0.0;
  setParameter("Pt_inlet",&Pt_inlet);
    
}

//////////////////////////////////////////////////////////////////////////////

SubOutletEuler2DPuvtMisMPI::~SubOutletEuler2DPuvtMisMPI()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler2DPuvtMisMPI::setup()
{
    
  FVMCC_BC::setup();
      
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler2DPuvtMisMPI::preProcess()
{
  const std::string nsp = getMethodData().getNamespace();
  _nP = PE::GetPE().GetProcessorCount(nsp);
  _myP = PE::GetPE().GetRank(nsp);
  _comm = PE::GetPE().GetCommunicator(nsp);
  
  //find the average pressure on the cells along the outlet boundary
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> _faceBuilder;
  _faceBuilder.setup();
  _faceBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  FaceTrsGeoBuilder::GeoData& faceData = _faceBuilder.getDataGE();
  faceData.isBFace = true;
  faceData.trs = trs;
  unsigned int nbTrsFaces = trs->getLocalNbGeoEnts();
  SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  CFuint n_Iter = subSysStatus->getNbIter();
  double P_sum = 0.0; 
  unsigned int n_adds = 0; 
  for (CFuint iFace = 0; iFace < nbTrsFaces; iFace++){
    CFLogDebugMed( "iFace = " << iFace << "\n");
    faceData.idx = iFace;
    GeometricEntity *const iface = _faceBuilder.buildGE();
    const CFuint faceGlobalID = iface->getID();
    State *const InnerState = iface->getState(0);
    if(InnerState->isParUpdatable()){
      P_sum = P_sum + (*InnerState)[0];
      n_adds = n_adds + 1;
    } 
    _faceBuilder.releaseGE();
  }
  _P_sum = P_sum;
  _Tot_P_sum = 0; 
  _n_adds = n_adds;
  _Tot_N_adds = _n_adds;
  MPI_Allreduce(&_P_sum, &_Tot_P_sum, 1, MPI_DOUBLE, MPI_SUM, _comm);
  MPI_Allreduce(&_n_adds, &_Tot_N_adds, 1, MPI_UNSIGNED, MPI_SUM, _comm);
  const CFreal gamma = (getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>())->getModel()->getGamma();
  _P_is = Pt_inlet/(pow((1 + ((gamma - 1)/2)*M_is*M_is),(gamma/(gamma - 1))));
  _P_average_inner = _Tot_P_sum/_Tot_N_adds;
  CFreal ratio_ave = _P_average_inner/_P_is;
  if(ratio_ave>1){
    ratio_ave = 1 / ratio_ave;
  }
  const CFreal damping_factor = 2;
  _damping = abs(pow((1 - abs(ratio_ave - 1)),damping_factor));
  std::cout <<"  P_is=   "<<_P_is  <<" P_average_inner=   "<<_P_average_inner  <<"  ratio_ave=   " <<ratio_ave  <<"  damping=   " <<_damping <<std::endl;
  
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletEuler2DPuvtMisMPI::setGhostState(GeometricEntity *const face)
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
  const CFreal M_inner = Vel/aInnerState;

  //subsonic outlet face case 
  if(M_normal_isoentropic<1.0){

  //std::cout <<"  P_is=   "<<_P_is  <<" P_average_inner=   "<<_P_average_inner  <<"  damping=   " <<_damping <<std::endl;    

    const CFreal P_exit = _P_is + ((*innerState)[0] - _P_average_inner)*_damping;
    const CFreal P_ghost = 2*P_exit - (*innerState)[0];
    
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

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
