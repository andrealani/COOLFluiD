#include "LaxFriedNSvtFlux.hh"
#include "Environment/ObjectProvider.hh"

#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/FVMCC_PolyRec.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<LaxFriedNSvtFlux,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeNavierStokesModule>
laxFriedNSvtFluxProvider("LaxFriedNSvt");

//////////////////////////////////////////////////////////////////////////////

LaxFriedNSvtFlux::LaxFriedNSvtFlux(const std::string& name) :
  LaxFriedFlux(name),
  m_nsVarSet(CFNULL),
  m_dummyGradients(),
  m_avState()
{
  addConfigOptionsTo(this);
  
  m_reynoldsMin = 0.1;
  setParameter("ReynoldsMin",&m_reynoldsMin);
  
  m_velocityIDs = vector<CFuint>();
  setParameter("VelocityIDs",&m_velocityIDs);
}

//////////////////////////////////////////////////////////////////////////////

LaxFriedNSvtFlux::~LaxFriedNSvtFlux()
{
}

//////////////////////////////////////////////////////////////////////////////

void LaxFriedNSvtFlux::setup()
{
  LaxFriedFlux::setup();
  
  m_nsVarSet = getMethodData().getDiffusiveVar().d_castTo
    <Physics::NavierStokes::NavierStokesVarSet>();
  
  m_avState.resize(PhysicalModelStack::getActive()->getNbEq());
  
  if (m_velocityIDs.size() == 0) {
    const CFuint dim = PhysicalModelStack::getActive()->getDim();
    m_velocityIDs.resize(dim);
    for (CFuint i = 0 ; i < dim; ++i) {
      m_velocityIDs[i] = 1 + i;
    }
  }
}
    
//////////////////////////////////////////////////////////////////////////////

CFreal LaxFriedNSvtFlux::getReductionCoeff()
{
  CellCenterFVMData& data = this->getMethodData(); 
  GeometricEntity& face = *data.getCurrentFace();
  const CFuint nbFaceNodes = face.nbNodes(); 
  DataHandle< RealVector> nstates = socket_nstates.getDataHandle();
    
  // average state based on extrapolated nodal values 
  m_avState = 0.0;
  for(CFuint i = 0; i < nbFaceNodes; ++i) {
     const CFuint nodeID = face.getNode(i)->getLocalID();
     m_avState += nstates[nodeID];
   }
  m_avState /= nbFaceNodes;
  
 // m_avState = 0.5*(data.getCurrRightState() + data.getCurrLeftState());
  
  const CFreal lstate = MathFunctions::getDistance
    (face.getState(LEFT)->getCoordinates(),
     face.getState(RIGHT)->getCoordinates());
  
  // ONLY 2D
  const CFreal lnode = MathFunctions::getDistance
   (*face.getNode(0), *face.getNode(1));
  
  const CFreal length = sqrt(lstate*lnode);
    
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  cf_assert(m_velocityIDs.size() == dim);
  
  CFreal speed = 0.0;
  for (CFuint i = 0; i < dim; ++i) {
    const CFreal avV = m_avState[m_velocityIDs[i]];
    speed += avV*avV;
  }
  speed = sqrt(speed);
  
  const CFreal nu = m_nsVarSet->getDynViscosity(m_avState, m_dummyGradients)/
    m_nsVarSet->getDensity(m_avState);
  
  //  cout << "speed  = " << speed << endl;
  //   cout << "length = " << length << endl;
  //   cout << "nu     = " << nu << endl;
  //if (speed*length/nu < 1.0) {	  
  // cout << "Reh    = " << speed*length/nu << endl << endl;
 // }

  // return the minimum between the user defined diffusion reduction 
  // coefficient and the local Reynolds number
  
  return min(_currentDiffRedCoeff, max(m_reynoldsMin,speed*length/nu));
}

//////////////////////////////////////////////////////////////////////////////

void LaxFriedNSvtFlux::defineConfigOptions(Config::OptionList& options)
{ 
  options.addConfigOption< CFreal >
    ("ReynoldsMin", "Minimum value for cell Reynolds number");
  
  options.addConfigOption< vector<CFuint> >
    ("VelocityIDs", "Variable IDs of the velocity components");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
