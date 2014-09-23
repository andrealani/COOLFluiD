#include "Framework/MethodStrategyProvider.hh"
#include "Framework/DataHandle.hh"
#include "Framework/GeometricEntity.hh"
#include "MHD/MHD3DProjectionVarSet.hh"
#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "FiniteVolumeMHD/MHD3DProjectionJupiterSourceTerm.hh"
#include "FiniteVolume/CellCenterFVMData.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<MHD3DProjectionJupiterSourceTerm,
               CellCenterFVMData,
               Framework::ComputeSourceTerm<CellCenterFVMData>,
               FiniteVolumeMHDModule>
mHD3DProjectionJupiterSTFVMCCProvider("MHD3DProjectionJupiterST");

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionJupiterSourceTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("distanceTorusCenterFromOrigin","Distance of center of the torus from the origin.");
  options.addConfigOption< CFreal >("radiusTorus","Radius of torus.");
  options.addConfigOption< CFreal >("velTangentialTorus","Corotation velocity of the torus (tangential).");
  options.addConfigOption< CFreal >("dRhodtTorusPlasma","mass density production rate in the torus.");
  options.addConfigOption< CFreal >("torusTemp","Temperature of the neutrals in the torus.");
  options.addConfigOption< CFreal >("torusAtomicMass","atomic mass of the ions in the torus.");
  options.addConfigOption< CFreal >("ionoRotPeriod","Rotation period of the planet.");
  options.addConfigOption< CFreal >("ionoCollFreq","Maximum Ion-Neutral collision frequency in the ionosphere.");
  options.addConfigOption< CFreal >("ionoDecay","Decay constant (it decreases exponentially) for the ionosphere.");
  options.addConfigOption< CFreal >("ionoAtomicMass","atomic mass of the ions in the ionosphere.");
  options.addConfigOption< CFreal >("radiusBnd","Radius of the inner boundary.");
  options.addConfigOption< CFreal >("ionoIonNeutralMassRatio","Mi/Mn - Mass ratio in the ionosphere.");
  options.addConfigOption< CFreal >("ionoTemp","Temperature of the neutrals in the ionosphere.");
  options.addConfigOption< CFreal >("gravityTimesJupiMass","G.M - Gravitational constant times the mass of the planet.");
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionJupiterSourceTerm::MHD3DProjectionJupiterSourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  _varSet(CFNULL),
  _physicalData(CFNULL),
  _stateCoord(CFNULL)
{
  addConfigOptionsTo(this);

  _a0 = 5.9;
  setParameter("distanceTorusCenterFromOrigin",&_a0);

  _radiusTorus = 1.0;
  setParameter("radiusTorus",&_radiusTorus);

  _velTangentialTorus = 1.0;
  setParameter("velTangentialTorus",&_velTangentialTorus);

  _dRhodtTorusPlasma = 1.0;
  setParameter("dRhodtTorusPlasma",&_dRhodtTorusPlasma);

  _torusTemp = 0.001;
  setParameter("torusTemp",&_torusTemp);

  _torusAtomicMass = 22.0;
  setParameter("torusAtomicMass",&_torusAtomicMass);
   
  _ionoRotPeriod = -149.92138;
  setParameter("ionoRotPeriod",&_ionoRotPeriod);

  _ionoCollFreq = 4885.15;
  setParameter("ionoCollFreq",&_ionoCollFreq);
    
  _ionoDecay = 0.25;
  setParameter("ionoDecay",&_ionoDecay);

  _ionoAtomicMass = 1.0;
  setParameter("ionoAtomicMass",&_ionoAtomicMass);
  
  _radiusBnd = 4.5;
  setParameter("radiusBnd",&_radiusBnd);
  	  
  _ionoIonNeutralMassRatio = 10.0;
  setParameter("ionoIonNeutralMassRatio",&_ionoIonNeutralMassRatio);
  	  
  _ionoTemp = 0.001;
  setParameter("ionoTemp",&_ionoTemp);
	  
  _gravityTimesJupiMass = 0.2228394;
  setParameter("gravityTimesJupiMass",&_gravityTimesJupiMass);
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionJupiterSourceTerm::~MHD3DProjectionJupiterSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionJupiterSourceTerm::setup()
{
  ComputeSourceTermFVMCC::setup();

  assert(_varSet.isNotNull());

  _varSet->getModel()->resizePhysicalData(_physicalData);
  
  _stateCoord.resize(Framework::PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionJupiterSourceTerm::setVarSet(Common::SafePtr<ConvectiveVarSet> varSet)
{
  _varSet = varSet.d_castTo<MHD3DProjectionVarSet>();
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionJupiterSourceTerm::computeSource(GeometricEntity *const element,
           RealVector& source)
{
  assert(_varSet.isNotNull());

  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  
  // this source term is for MHD flows
  const vector<State*>* const states = element->getStates();
  const CFuint elementID = element->getID();

  // all elements in FVM should have only one state
  assert(states->size() == 1);

  State *const currState = (*states)[0];

  _stateCoord = currState->getCoordinates();
  const CFreal mX = _varSet->getModel()->getMX();
  const CFreal mY = _varSet->getModel()->getMY();
  const CFreal mZ = _varSet->getModel()->getMZ();

  _varSet->computePhysicalData(*currState,_physicalData);
  const CFreal pressure = _physicalData[MHDTerm::P];

  // distance from the centre of the domain
  const CFreal radius = sqrt(_stateCoord[0]*_stateCoord[0]+_stateCoord[1]*_stateCoord[1]+_stateCoord[2]*_stateCoord[2]);
 
  const CFreal pi = 3.14159265358979323846;
  const CFreal twoPi = 2.0*pi;
  
  // velocity of the neutrals in the ionosphere
  CFreal vnX=0.0, vnY=0.0, vnZ=0.0;
  
  // velocity of the neutrals in the ionosphere 
  // for all the possible orientations of the dipole moment
  if ((mX == 0.0) && (mY == 0.0)) {
    const CFreal twoPi_T=twoPi/_ionoRotPeriod*mZ/fabs(mZ);
    vnX = twoPi_T*_stateCoord[1];
    vnY = -twoPi_T*_stateCoord[0];
    vnZ = 0.0;
    }
  else if ((mX > 0.0) && (mY == 0.0)) {
    const CFreal theta = atan(mZ/mX);
    const CFreal twoPi_T = twoPi/_ionoRotPeriod;
    vnX = twoPi_T*_stateCoord[1]*sin(theta);
    vnY = twoPi_T*(_stateCoord[2]*cos(theta)-_stateCoord[0]*sin(theta));
    vnZ = -twoPi_T*_stateCoord[1]*cos(theta);
    }
  else if ((mX < 0.0) && (mY == 0.0)) {
    const CFreal theta = pi + atan(mZ/mX);
    const CFreal twoPi_T = twoPi/_ionoRotPeriod;
    vnX = twoPi_T*_stateCoord[1]*sin(theta);
    vnY = twoPi_T*(_stateCoord[2]*cos(theta)-_stateCoord[0]*sin(theta));
    vnZ = -twoPi_T*_stateCoord[1]*cos(theta);
    }
  else if (mY > 0.0) {
    const CFreal theta = atan(mZ/sqrt(mX*mX+mY*mY));
    const CFreal phi = atan(mX/mY);
    const CFreal twoPi_T = twoPi/_ionoRotPeriod;
    vnX = twoPi_T*(_stateCoord[1]*sin(theta)-_stateCoord[2]*cos(theta)*cos(phi));
    vnY = twoPi_T*(_stateCoord[2]*cos(theta)*sin(phi)-_stateCoord[0]*sin(theta));
    vnZ = -twoPi_T*(_stateCoord[0]*cos(theta)*cos(phi)-_stateCoord[1]*sin(phi)*cos(theta));
    }
  else if (mY < 0.0) {
    const CFreal theta = atan(mZ/sqrt(mX*mX+mY*mY));
    const CFreal phi = pi + atan(mX/mY);
    const CFreal twoPi_T = twoPi/_ionoRotPeriod;
    vnX = twoPi_T*(_stateCoord[1]*sin(theta)-_stateCoord[2]*cos(theta)*cos(phi));
    vnY = twoPi_T*(_stateCoord[2]*cos(theta)*sin(phi)-_stateCoord[0]*sin(theta));
    vnZ = -twoPi_T*(_stateCoord[0]*cos(theta)*cos(phi)-_stateCoord[1]*sin(phi)*cos(theta));
    }
    
  // the collision frequency decrease exponentially with the radial distance
  // because the density of the neutrals decrease exponentially with the radial distance
  const CFreal nu_in = _ionoCollFreq*exp((_radiusBnd-radius)/_ionoDecay);
  const CFreal rhoNu_in = (*currState)[0]*nu_in;
  
  const CFreal vx = (*currState)[1]/(*currState)[0];
  const CFreal vy = (*currState)[2]/(*currState)[0];
  const CFreal vz = (*currState)[3]/(*currState)[0];

  // kinetic energy due to the collisions between the neutrals and the ions in the ionosphere
  const CFreal energy1 = 1.0/(1.0+_ionoIonNeutralMassRatio)*nu_in*((vnX-vx)*(vnX-vx)+(vnY-vy)*(vnY-vy)+(vnZ-vz)*(vnZ-vz))*(*currState)[0];

  // thermal energy due to the collisions between the neutrals and the ions in the ionosphere	
  // 8250.6277 is k/m_h
  // k/m_i=k/(A*m_h) 
  const CFreal energy2 = (8250.6277/_ionoAtomicMass*_ionoTemp*(*currState)[0]-pressure)*3.0/(1.0+1.0/_ionoIonNeutralMassRatio)*nu_in;
  
  // extra thermal energy due to the mass loading in the torus
  CFreal energy3 = 0.0;
  
  // velocity of the neutrals in the torus (uTorusPlasma,vTorusPlasma,wTorusPlasma)	  
  // mass loading rate (rhoChangeTorusPlasma)
  // extra kinetic energy due to the mass loading (extra)
  CFreal uTorusPlasma = 0.0, vTorusPlasma = 0.0, wTorusPlasma = 0.0, rhoChangeTorusPlasma = 0.0, extra=0.0;
  
  // calculation of the extra thermal and kinetic energy 
  // and of the velocity of the neutrals in the torus
  // for all the possible orientations of the dipole moment
  if ((mX == 0.0) && (mY == 0.0)) {
    const CFreal a = sqrt(_stateCoord[0]*_stateCoord[0]+_stateCoord[1]*_stateCoord[1]);
    const CFreal distTorusCentre = sqrt((a-_a0)*(a-_a0) + _stateCoord[2]*_stateCoord[2]);
    if (distTorusCentre <= _radiusTorus) {
      const CFreal vt_a = _velTangentialTorus/a*mZ/fabs(mZ);
      uTorusPlasma = vt_a*_stateCoord[1];
      vTorusPlasma = -vt_a*_stateCoord[0];
      rhoChangeTorusPlasma = _dRhodtTorusPlasma*0.5*(1.0-cos(pi*(_radiusTorus-distTorusCentre)/_radiusTorus));
      extra=0.5*rhoChangeTorusPlasma*(uTorusPlasma*uTorusPlasma+vTorusPlasma*vTorusPlasma+wTorusPlasma*wTorusPlasma);
      // 8250.6277 is k/m_h
      //  k/m_i=k/(A*m_h)
      //  8250.6277 *1.5 = 12375.942
      energy3 =12375.942/_torusAtomicMass*_torusTemp*rhoChangeTorusPlasma;
      }
  }
  else if ((mZ == 0.0) && (mY == 0.0) && (mX != 0.0)) {
    const CFreal a = sqrt(_stateCoord[2]*_stateCoord[2]+_stateCoord[1]*_stateCoord[1]);
    const CFreal distTorusCentre = sqrt((a-_a0)*(a-_a0) + _stateCoord[0]*_stateCoord[0]);
    if (distTorusCentre <= _radiusTorus) {
      const CFreal vt_a = _velTangentialTorus/a*mZ/fabs(mZ);
      vTorusPlasma = vt_a*_stateCoord[2];
      wTorusPlasma = -vt_a*_stateCoord[1];
      rhoChangeTorusPlasma = _dRhodtTorusPlasma*0.5*(1.0-cos(pi*(_radiusTorus-distTorusCentre)/_radiusTorus));
      extra=0.5*rhoChangeTorusPlasma*(uTorusPlasma*uTorusPlasma+vTorusPlasma*vTorusPlasma+wTorusPlasma*wTorusPlasma);
      // 8250.6277 is k/m_h
      // k/m_i=k/(A*m_h)
      // 8250.6277 *1.5 = 12375.942
      energy3 =12375.942/_torusAtomicMass*_torusTemp*rhoChangeTorusPlasma;
      }
    }
  else if ((mY > 0.0)) {
    const CFreal theta = atan(mZ/sqrt(mX*mX+mY*mY));
    const CFreal phi = atan(mX/mY) ;
    const CFreal newX = _stateCoord[0]*cos(phi) - _stateCoord[1]*sin(phi);
    const CFreal newY = _stateCoord[1]*cos(phi)*sin(theta) + _stateCoord[0]*sin(phi)*sin(theta) - _stateCoord[2]*cos(theta);
    const CFreal newZ = _stateCoord[1]*cos(phi)*cos(theta) + _stateCoord[0]*sin(phi)*cos(theta) + _stateCoord[2]*sin(theta);
    const CFreal a = sqrt(newX*newX+newY*newY);
    const CFreal distTorusCentre = sqrt((a-_a0)*(a-_a0) + newZ*newZ);
    if (distTorusCentre <= _radiusTorus) {
      const CFreal vt_a = _velTangentialTorus/a;
      const CFreal vNewX = vt_a*newY;
      const CFreal vNewY = -vt_a*newX;
      uTorusPlasma = vNewX*cos(phi) + vNewY*sin(theta)*sin(phi);
      vTorusPlasma = vNewY*sin(phi)*cos(theta) - vNewX*sin(phi);
      wTorusPlasma = -vNewY*cos(theta);
      rhoChangeTorusPlasma = _dRhodtTorusPlasma*0.5*(1.0-cos(pi*(_radiusTorus-distTorusCentre)/_radiusTorus));
      extra=0.5*rhoChangeTorusPlasma*(uTorusPlasma*uTorusPlasma+vTorusPlasma*vTorusPlasma+wTorusPlasma*wTorusPlasma);
      // 8250.6277 is k/m_h
      // k/m_i=k/(A*m_h)
      // 8250.6277 *1.5 = 12375.942
      energy3 =12375.942/_torusAtomicMass*_torusTemp*rhoChangeTorusPlasma;
      }
    }
  else if ((mY < 0.0) && (mX != 0.0)) {
    const CFreal theta = atan(mZ/sqrt(mX*mX+mY*mY));
    const CFreal phi = atan(mX/mY) + pi*fabs(mX)/mX;
    const CFreal newX = _stateCoord[0]*cos(phi) - _stateCoord[1]*sin(phi);
    const CFreal newY = _stateCoord[1]*cos(phi)*sin(theta) + _stateCoord[0]*sin(phi)*sin(theta) - _stateCoord[2]*cos(theta);
    const CFreal newZ = _stateCoord[1]*cos(phi)*cos(theta) + _stateCoord[0]*sin(phi)*cos(theta) + _stateCoord[2]*sin(theta);
    const CFreal a = sqrt(newX*newX+newY*newY);
    const CFreal distTorusCentre = sqrt((a-_a0)*(a-_a0) + newZ*newZ);
    if (distTorusCentre <= _radiusTorus) {
      const CFreal vt_a = _velTangentialTorus/a;
      const CFreal vNewX = vt_a*newY;
      const CFreal vNewY = -vt_a*newX;
      uTorusPlasma = vNewX*cos(phi) + vNewY*sin(theta)*sin(phi);
      vTorusPlasma = vNewY*sin(phi)*cos(theta) - vNewX*sin(phi);
      wTorusPlasma = -vNewY*cos(theta);
      rhoChangeTorusPlasma = _dRhodtTorusPlasma*0.5*(1.0-cos(pi*(_radiusTorus-distTorusCentre)/_radiusTorus));
      extra=0.5*rhoChangeTorusPlasma*(uTorusPlasma*uTorusPlasma+vTorusPlasma*vTorusPlasma+wTorusPlasma*wTorusPlasma);
      // 8250.6277 is k/m_h
      // k/m_i=k/(A*m_h)
      // 8250.6277 *1.5 = 12375.942
      energy3 =12375.942/_torusAtomicMass*_torusTemp*rhoChangeTorusPlasma;
      }
    }
  else if ((mY < 0.0) && (mX == 0.0)) {
    const CFreal theta = atan(mZ/mY);
    const CFreal newX = -_stateCoord[0];
    const CFreal newY = _stateCoord[1]*sin(theta) - _stateCoord[2]*cos(theta);
    const CFreal newZ = _stateCoord[1]*cos(theta) + _stateCoord[2]*sin(theta);
    const CFreal a = sqrt(newX*newX+newY*newY);
    const CFreal distTorusCentre = sqrt((a-_a0)*(a-_a0) + newZ*newZ);
    if (distTorusCentre <= _radiusTorus) {
      const CFreal vt_a = _velTangentialTorus/a;
      uTorusPlasma = -vt_a*newY;
      vTorusPlasma = -sin(theta)*vt_a*newX;
      wTorusPlasma = cos(theta)*vt_a*newX;
      rhoChangeTorusPlasma = _dRhodtTorusPlasma*0.5*(1.0-cos(pi*(_radiusTorus-distTorusCentre)/_radiusTorus));
      extra=0.5*rhoChangeTorusPlasma*(uTorusPlasma*uTorusPlasma+vTorusPlasma*vTorusPlasma+wTorusPlasma*wTorusPlasma);
      // 8250.6277 is k/m_h
      // k/m_i=k/(A*m_h)
      // 8250.6277 *1.5 = 12375.942
      energy3 =12375.942/_torusAtomicMass*_torusTemp*rhoChangeTorusPlasma;
      }
    }
  else if ((mY == 0.0) && (mX != 0.0) && (mZ != 0.0)) {
    const CFreal theta = fabs(mX)/mX*1.570796327 - atan(mZ/mY);
    const CFreal newX = _stateCoord[0]*cos(theta) - _stateCoord[2]*sin(theta);
    const CFreal newY = _stateCoord[1];
    const CFreal newZ = _stateCoord[0]*sin(theta) + _stateCoord[2]*cos(theta);
    const CFreal a = sqrt(newX*newX+newY*newY);
    const CFreal distTorusCentre = sqrt((a-_a0)*(a-_a0) + newZ*newZ);
    if (distTorusCentre <= _radiusTorus) {
      const CFreal vt_a = _velTangentialTorus/a;
      uTorusPlasma = cos(theta)*vt_a*newY;
      vTorusPlasma = -vt_a*newX;
      wTorusPlasma = -sin(theta)*vt_a*newY;
      rhoChangeTorusPlasma = _dRhodtTorusPlasma*0.5*(1.0-cos(pi*(_radiusTorus-distTorusCentre)/_radiusTorus));
      extra=0.5*rhoChangeTorusPlasma*(uTorusPlasma*uTorusPlasma+vTorusPlasma*vTorusPlasma+wTorusPlasma*wTorusPlasma);
      // 8250.6277 is k/m_h
      // k/m_i=k/(A*m_h)
      // 8250.6277 *1.5 = 12375.942
      energy3 =12375.942/_torusAtomicMass*_torusTemp*rhoChangeTorusPlasma;
      }
    }
	  
  const CFreal mMG_Rcube = -_gravityTimesJupiMass/(radius*radius*radius);
  const CFreal mMGRho_Rcube = (*currState)[0]*mMG_Rcube;

  source[0] = rhoChangeTorusPlasma*volumes[elementID];
  source[1] = (rhoNu_in*(vnX-vx)+rhoChangeTorusPlasma*uTorusPlasma+mMGRho_Rcube*_stateCoord[0])*volumes[elementID];
  source[2] = (rhoNu_in*(vnY-vy)+rhoChangeTorusPlasma*vTorusPlasma+mMGRho_Rcube*_stateCoord[1])*volumes[elementID];
  source[3] = (rhoNu_in*(vnZ-vz)+rhoChangeTorusPlasma*wTorusPlasma+mMGRho_Rcube*_stateCoord[2])*volumes[elementID];
  source[4] = 0.0;
  source[5] = 0.0;
  source[6] = 0.0;
  source[7] = (extra+mMG_Rcube*((*currState)[1]*_stateCoord[0]+(*currState)[2]*_stateCoord[1]+(*currState)[3]*_stateCoord[2])+energy1+energy2+energy3)*volumes[elementID];
  source[8] = 0.0;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
