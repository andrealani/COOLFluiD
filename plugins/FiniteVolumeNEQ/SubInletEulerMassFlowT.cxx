#include "MathTools/MathConsts.hh"
#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "FiniteVolumeNEQ/SubInletEulerMassFlowT.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MultiScalarTerm.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Common/BadValueException.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SubInletEulerMassFlowT, CellCenterFVMData, FiniteVolumeNEQModule>
subInletEuler2DMassFlowTFVMCCProvider("SubInletEulerMassFlowTFVMCC");

//////////////////////////////////////////////////////////////////////////////

SubInletEulerMassFlowT::SubInletEulerMassFlowT(const std::string& name) :
  SubInletEulerFunc(name),
  m_term(CFNULL),
  m_library(CFNULL),
  m_stateHasPartialDensities(true),
  m_masses()
{
  addConfigOptionsTo(this);
  
  m_yi = std::vector<CFreal>();
  setParameter("Yi",&m_yi);
}

//////////////////////////////////////////////////////////////////////////////

SubInletEulerMassFlowT::~SubInletEulerMassFlowT()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubInletEulerMassFlowT::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<CFreal> >("Yi","Inlet mass fractions.");
}
      
//////////////////////////////////////////////////////////////////////////////

void SubInletEulerMassFlowT::setGhostState(GeometricEntity *const face)
{
  CFLog(DEBUG_MIN, "SubInletEulerMassFlowT::setGhostState()\n");

  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);	
  
  if(_useFunction)
  {
    // coordinate of the boundary point
    _bCoord = (innerState->getCoordinates() + ghostState->getCoordinates());
    _bCoord *= 0.5;
    
    // (*ghostState) = 2*bcState - (*innerState)
    _vFunction.evaluate(_bCoord, _inletData);
    
    _inletData[1] /= m_term->getTempRef();
  }
  else {
    _inletData[0] = _massFlow;
    _inletData[1] = _temperature;
  }
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbSpecies = m_library->getNbSpecies();
  const CFuint nbTs = m_library->getNbTempVib() + m_library->getNbTe();
  const CFuint TID = nbSpecies + dim;
  
  this->getMethodData().getUpdateVar()->computePhysicalData(*innerState, _dataInnerState);
  
  // note that we impose T=Te in the inlet, otherwise, if T!=Te, we would need to 
  // consider p=p_h(T)+p_e(T_e) which would make things a little different
  const CFreal Tinlet = _inletData[1];
  cf_assert(Tinlet > 0.);
  
  const CFreal Tin = 2.*Tinlet - (*innerState)[TID];
  (*ghostState)[TID] = (Tin > 0.) ? Tin : Tinlet; 
  cf_assert((*ghostState)[TID] > 0.);
  
  for (CFuint i = 0; i < nbTs; ++i) {
    const CFuint tvID = TID+1+i;
    const CFreal TinI = 2.*Tinlet - (*innerState)[tvID];
    (*ghostState)[tvID] = (TinI > 0.) ? TinI : Tinlet; 
    cf_assert((*ghostState)[tvID] > 0.);
  }
  
  CFreal ovMMass = 0.;
  for (CFuint i = 0; i < nbSpecies; ++i) {
    ovMMass += m_yi[i]/m_masses[i];
  }
  cf_assert(ovMMass > 0.);
  
  const CFreal RovM = m_library->getRgas()*ovMMass;
  // p_inner = p_inlet = rho_inlet * R * T_inlet / M_inlet(y_inlet) 
  const CFreal pinner = m_term->getPressureFromState(_dataInnerState[EulerTerm::P]);
  cf_assert(pinner > 0.);
  const CFreal rho  = pinner/(RovM*Tinlet);
  cf_assert(rho > 0.);
  const CFreal area = MathTools::MathConsts::CFrealPi()*(_inletRadii[1]*_inletRadii[1] - 
							 _inletRadii[0]*_inletRadii[0]);
  const CFreal uInf = .001*_inletData[0]/(area*rho); // mass flow is given in g/s instead of Kg/s
  cf_assert(uInf > 0.);
  
  // we extrapolate from inside the density (pressure)
  // we impose y_i, velocity and temperatures
  if (!m_stateHasPartialDensities) {
    // pressure is extrapolated from inside so we impose gradient=0 for partial pressures
    // for (CFuint i = 0; i < nbSpecies; ++i) {
    //  (*ghostState)[i] = (*innerState)[i];
    // }
    
    const CFreal rhoG = 2.*rho - _dataInnerState[EulerTerm::RHO];    
    CFreal sumY = 0.;
    const CFuint startSpecies = m_term->getFirstScalarVar(0);
    for (CFuint i = 0; i < nbSpecies; ++i) {
      const CFreal yi = 2.*m_yi[i] - _dataInnerState[startSpecies+i];
      sumY += yi; 
      // AL: check this better ...
      const CFreal Ri = m_library->getRgas()/m_masses[i];
     
      CFreal Ti = (*ghostState)[TID];  
      if (nbTs > 0) {  
	Ti = (i > 0) ? (*ghostState)[TID] : (*ghostState)[TID+nbTs]; // Te is always the last temperature
      }
      (*ghostState)[i] = rhoG*yi*Ri*Ti; // rho_i
    } 
    
    // sanity check: sum(y_i)=1
    cf_assert(sumY > 0.9999 && sumY < 1.0001);
  }
  else {
    const CFreal rhoG = 2.*rho - _dataInnerState[EulerTerm::RHO];    
    CFreal sumY = 0.;
    const CFuint startSpecies = m_term->getFirstScalarVar(0);
    for (CFuint i = 0; i < nbSpecies; ++i) {
      const CFreal yi = 2.*m_yi[i] - _dataInnerState[startSpecies+i];
      sumY += yi;
      (*ghostState)[i] = rhoG*yi; // rho_i
    } 
    
    // sanity check: sum(y_i)=1
    cf_assert(sumY > 0.9999 && sumY < 1.0001);
  }
  
  (*ghostState)[nbSpecies]   = 2.*uInf - _dataInnerState[EulerTerm::VX];
  (*ghostState)[nbSpecies+1] = - _dataInnerState[EulerTerm::VY]; // 0; here there was a huge bug
  if (dim == DIM_3D) {(*ghostState)[nbSpecies+2] = - _dataInnerState[EulerTerm::VZ];}
  
  CFLog(DEBUG_MAX, "SubInletEulerMassFlowT::setGhostState() => rho_i u v T Tv ER EI = " << *ghostState << "\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void SubInletEulerMassFlowT::setup()
{
  SubInletEulerFunc::setup();
  
  m_term = PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm().d_castTo<MultiScalarTerm<EulerTerm> >();
  
  m_library = PhysicalModelStack::getActive()->
    getImplementor()->getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  const vector<string>& varNames = this->getMethodData().getUpdateVar()->getVarNames();
  if ((int) std::count(varNames.begin(), varNames.end(), "rho0") > 0) {
    m_stateHasPartialDensities = true;
  }
  else if ((int) std::count(varNames.begin(), varNames.end(), "p0") > 0) {
    m_stateHasPartialDensities = false;
  }
  
  const CFuint nbSpecies = m_library->getNbSpecies();
  m_masses.resize(nbSpecies);
  m_library->getMolarMasses(m_masses);
  
  for (CFuint i = 0; i < nbSpecies; ++i) {
    cf_assert(m_yi[i] >= 0.);
  }
} 

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
