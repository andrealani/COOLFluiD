#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"

#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/Euler2DVarSet.hh"

#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/SubOutletEulerPuvtMisMPI.hh"

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

MethodCommandProvider<SubOutletEulerPuvtMisMPI, CellCenterFVMData, FiniteVolumeNavierStokesModule>
subOutletEulerPuvtMisMPIFVMCCProvider("SubOutletEulerPuvtMisMPIFVMCC");

//////////////////////////////////////////////////////////////////////////////

void SubOutletEulerPuvtMisMPI::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal > ("Mis","M_is");
  options.addConfigOption< CFreal > ("Pt_inlet","Pt_inlet");
}

//////////////////////////////////////////////////////////////////////////////

SubOutletEulerPuvtMisMPI::SubOutletEulerPuvtMisMPI(const std::string& name) :
  FVMCC_BC(name)
{
  addConfigOptionsTo(this);
  M_is = 0.0;
  setParameter("Mis",&M_is);
  Pt_inlet = 0.0;
  setParameter("Pt_inlet",&Pt_inlet);
    
}

//////////////////////////////////////////////////////////////////////////////

SubOutletEulerPuvtMisMPI::~SubOutletEulerPuvtMisMPI()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletEulerPuvtMisMPI::setup()
{
    
  FVMCC_BC::setup();
      
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletEulerPuvtMisMPI::preProcess()
{
  const std::string nsp = this->getMethodData().getNamespace();
  _nP = PE::GetPE().GetProcessorCount(nsp);
  _myP = PE::GetPE().GetRank(nsp);
  _comm = PE::GetPE().GetCommunicator(nsp);
  
  _nbDim = PhysicalModelStack::getActive()->getDim();
  //find the average pressure on the cells along the outlet boundary
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> _faceBuilder;
  _faceBuilder.setup();
  _faceBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  FaceTrsGeoBuilder::GeoData& faceData = _faceBuilder.getDataGE();
  faceData.isBFace = true;
  faceData.trs = trs;
  unsigned int nbTrsFaces = trs->getLocalNbGeoEnts();
  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  
  //_Pghost.reserve(nbTrsFaces,_P_is);
  double P_sum = 0.0; 
  double P_sum_ghost = 0.0; 
  unsigned int n_adds = 0; 
  for (CFuint iFace = 0; iFace < nbTrsFaces; iFace++){
    CFLogDebugMed( "iFace = " << iFace << "\n");
    faceData.idx = iFace;
    GeometricEntity *const iface = _faceBuilder.buildGE();
    State *const InnerState = iface->getState(0);
    State *const GhostState = iface->getState(1);
    if(InnerState->isParUpdatable()){
      P_sum = P_sum + (*InnerState)[0];
      P_sum_ghost = P_sum_ghost + (*GhostState)[0];
      n_adds = n_adds + 1;
      //cout<<" (*InnerState)[0] = "<<(*InnerState)[0]<<endl;
    } 
    _faceBuilder.releaseGE();
  }
  _P_sum = P_sum;
  _Tot_P_sum = 0; 
  CFreal _P_sum_ghost = P_sum_ghost;
  CFreal _Tot_P_sum_ghost = 0; 
  _n_adds = n_adds;
  _Tot_N_adds = _n_adds;
  MPI_Allreduce(&_P_sum, &_Tot_P_sum, 1, MPI_DOUBLE, MPI_SUM, _comm);
  MPI_Allreduce(&_P_sum_ghost, &_Tot_P_sum_ghost, 1, MPI_DOUBLE, MPI_SUM, _comm);
  MPI_Allreduce(&_n_adds, &_Tot_N_adds, 1, MPI_UNSIGNED, MPI_SUM, _comm);

  CFreal gamma = 0.;
  if(_nbDim == 2){
    gamma = (getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>())->getModel()->getGamma();
  }
  if(_nbDim == 3){
    gamma = (getMethodData().getUpdateVar().d_castTo<Euler3DVarSet>())->getModel()->getGamma();
  }
  _P_is = Pt_inlet/(pow((1 + ((gamma - 1)/2)*M_is*M_is),(gamma/(gamma - 1))));
  _P_average_inner = _Tot_P_sum/_Tot_N_adds;
  _P_average_ghost = _Tot_P_sum_ghost/_Tot_N_adds;
  CFreal ratio_ave = _P_average_inner/_P_is;

  /*CFreal pgg = _P_average_ghost/_P_average_inner;
  CFreal pratio = abs(1 - pgg);
  if(pratio > 1){
    pratio = 1/pratio;
  }
  if(relaxationFactor != 1){
  relaxationFactor = 0.1;
  if(pratio<0.0001){
    relaxationFactor = 1;
  }
  }*/
  if(ratio_ave>1){
    ratio_ave = 1 / ratio_ave;
  }
  const CFreal damping_factor = 1;
  _damping = abs(pow((1 - abs(ratio_ave - 1)),damping_factor));
  _globalToLocalTRSFaceID.reserve(nbTrsFaces);
  _ratio1.clear();
  _ratio2.clear();
  _ghost.clear();
  _ratio1.reserve(nbTrsFaces);
  _ratio2.reserve(nbTrsFaces);
  _ghost.reserve(nbTrsFaces);
  CFreal as = 0;
  for (CFuint iFace = 0; iFace < nbTrsFaces; iFace++){
    faceData.idx = iFace;
    GeometricEntity *const iface = _faceBuilder.buildGE();
    const CFuint faceGlobalID = iface->getID();
    State *const InnerState = iface->getState(0);
    State *const GhostState = iface->getState(1);
    _globalToLocalTRSFaceID.insert(faceGlobalID,iFace);
    CFreal Ratio1 = (*InnerState)[0]/(*GhostState)[0];
    CFreal Ratio2 = _P_average_inner/(*InnerState)[0];
    CFreal a = Ratio1*(2*_damping - 1) - 2*_damping*Ratio1*Ratio2;
    if(a > as){
      as = a;
    }
    CFreal pgg = (*GhostState)[0]/_ghost[iFace];
    CFreal pratio = abs(1 - pgg);
    if(pratio > 1){
      pratio = 1/pratio;
    }
    _ratio1.push_back(Ratio1);
    _ratio2.push_back(Ratio2);
    _ghost.push_back((*GhostState)[0]);
    _faceBuilder.releaseGE();
  }
  asn2 = 1;
  if(abs(as)>=1){
    asn2 = 1/(abs(as) + 0.1);
  }  
  std::cout <<"  P_is=   "<<_P_is  <<" P_average_inner=   "<<_P_average_inner  <<" _P_average_ghost = "<<_P_average_ghost<<"  ratio_ave=   " <<ratio_ave  <<"  damping=   " <<_damping <<std::endl;
  
}

//////////////////////////////////////////////////////////////////////////////

void SubOutletEulerPuvtMisMPI::setGhostState(GeometricEntity *const face)
{
    
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  CFreal gamma = 0.;
  CFreal R = 0.;
  if(_nbDim == 2){
    gamma = (getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>())->getModel()->getGamma();
    R = (getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>())->getModel()->getR();
  }
  if(_nbDim == 3){
   gamma = (getMethodData().getUpdateVar().d_castTo<Euler3DVarSet>())->getModel()->getGamma();
   R = (getMethodData().getUpdateVar().d_castTo<Euler3DVarSet>())->getModel()->getR();
  }
  CFreal aInnerState = 0.;
  CFreal u = (*innerState)[1];
  CFreal v = (*innerState)[2];
  CFreal w = 0.;
  const CFuint faceID = face->getID();
  const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  CFreal nx = normals[startID];
  CFreal ny = normals[startID + 1];
  CFreal nz = 0.;
  if(_nbDim == 2){
    nz = 0;
    w = 0;
    aInnerState = sqrt(gamma*R*(*innerState)[3]);
  }
  if(_nbDim == 3){
    nz = normals[startID + 2];
    w = (*innerState)[3];
    aInnerState = sqrt(gamma*R*(*innerState)[4]);
  }
  CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny + nz*nz);
  nx *= invFaceLength;
  ny *= invFaceLength;
  nz *= invFaceLength;
  
  const CFreal vn = u*nx + v*ny + w*nz;
  //  const CFreal Vel = sqrt(u*u + v*v + w*w);
  const CFreal M_normal_isoentropic = vn / aInnerState;
  // const CFreal M_inner = Vel/aInnerState;

  //subsonic outlet face case 
  if(M_normal_isoentropic<1.0){

    CFreal iFace = _globalToLocalTRSFaceID.find(faceID);
    CFreal a = _ratio1[iFace]*(2*_damping - 1) - 2*_damping*_ratio1[iFace]*_ratio2[iFace];
    CFreal asn;
    if(abs(a)<1){
      asn = 1;
    }
    if(abs(a)>=1){
      asn = 1/(abs(a) + 0.1);
    }
    //cout<<" a = "<<a<<" asn = "<<asn<<" asn2 = "<<asn2<<endl;
    CFreal ConvFactCoeff2 = 0.8;
    if((*innerState)[0] > _P_average_inner){
      ConvFactCoeff2 = 1.2;
    }
    CFreal ConvFactCoeff = 0;
    if(MathTools::MathChecks::isEqualWithError((*innerState)[0],_P_average_inner, 0.01)){
      ConvFactCoeff = 0.008;
      ConvFactCoeff2 = 0.0;
    }
    //    CFreal ConvFact = ConvFactCoeff2*(((*innerState)[0] - 2*_P_is)/(ConvFactCoeff + 2*((*innerState)[0] - _P_average_inner)));
    CFreal P_exit = _P_is + ((*innerState)[0] - _P_average_inner)*_damping;
    //CFreal P_exit = (*innerState)[0]*(_P_is/_P_average_inner);
    //CFreal P_exit = _P_is + (*innerState)[0] - _P_average_ghost;
    CFreal P_ghost = 2*P_exit - (*innerState)[0];
    // CFreal isG = abs((_ghost[iFace] - P_ghost)/_ghost[iFace]);
    //if(MathTools::MathChecks::isGreaterWithError(isG, 0.25, 0.05)){
      //P_ghost = (_P_is + _ghost[iFace])/2;
    //}
    //cout<<" ConvFact = "<<ConvFact<<" (*innerState)[0] = "<<(*innerState)[0]<<" P_exit = "<<P_exit<<" P_ghost = "<<P_ghost<<endl;
    /*CFreal t = 0.1;
    CFreal fuzz = (abs(P_ghost - (*innerState)[0]))/abs(P_ghost);
    if(fuzz >=t){
      if(P_ghost > (*innerState)[0]){
        P_ghost = (*innerState)[0]/(1 - t);
      }
      if(P_ghost < (*innerState)[0]){
        P_ghost = (*innerState)[0]/(1 + t);
      }
    }
    if(P_ghost > Pt_inlet){
      P_ghost = Pt_inlet;
    }
    if(P_ghost < 50000){
      P_ghost = 50000;
    }
    /*
    CFreal pig = P_ghost/(*innerState)[0];
    CFreal pgg = P_ghost/_ghost[iFace];
    CFreal pratio = abs(1 - pgg);
    if(pratio > 1){
      pratio = 1/pratio;
    }
    //relaxationFactor = exp(-100*abs(pratio));
    /*if(relaxationFactor<0.1){
      relaxationFactor = 0.1;
    }*/
    /*relaxationFactor = 1.5;
    /*if(MathTools::MathChecks::isEqualWithError(pgg,1, 0.001)){
      relaxationFactor = 1;
    }
    /*if(pratio<0.0005){
      relaxationFactor = 1;
    }*/
    //cout<<" ratio g/i = "<<pig<<" relaxationFactor = "<<relaxationFactor<<endl;
    /*P_ghost  = relaxationFactor*P_ghost + (1 - relaxationFactor)*_ghost[iFace];*/
    //cout<<" (*innerState)[0] "<<(*innerState)[0]<<" P_exit = "<<P_exit<<" P_ghost = "<<P_ghost<<endl;
    //if(_myP==0){
      //P_ghost = _P_is;
      //(*innerState)[0] = _P_is;
      //cout<<" WARNING GHOST PRESSURE "<<endl;
      //cout<<" faceID "<<faceID<<" P_ghost = "<<P_ghost<<" (*innerState)[0] = "<<(*innerState)[0]<<endl;

    //}
    
    // set the variables in the ghost cell
    (*ghostState)[0] = P_ghost;
    (*ghostState)[1] = (*innerState)[1];
    (*ghostState)[2] = (*innerState)[2];
    (*ghostState)[3] = (*innerState)[3];
    if(_nbDim == 3){
     (*ghostState)[4] = (*innerState)[4];
    }
/*
    if((*innerState)[0] > Pt_inlet){
      cout<<" WARNING INNER PRESSURE SUBSONIC"<<endl;
    }*/

  }
  

  //supersonic outlet face case
  if(M_normal_isoentropic>=1.0){
    CFreal P_ghost2 = (*innerState)[0];
    if(P_ghost2 > Pt_inlet){
      P_ghost2 = Pt_inlet;
    }
    if(P_ghost2 < 50000){
      P_ghost2 = 50000;
    }
    
    (*ghostState)[0] = (*innerState)[0];    
    (*ghostState)[1] = (*innerState)[1];
    (*ghostState)[2] = (*innerState)[2];
    (*ghostState)[3] = (*innerState)[3];
    if(_nbDim == 3){
     (*ghostState)[4] = (*innerState)[4];
    }

  }

}


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
