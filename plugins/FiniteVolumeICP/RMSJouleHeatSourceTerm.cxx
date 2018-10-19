#include "FiniteVolumeICP/RMSJouleHeatSourceTerm.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "NavierStokes/EulerTerm.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "Framework/PhysicalConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<RMSJouleHeatSourceTerm,
		       CellCenterFVMData, 
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeModule>
rmsJouleHeatSTFVMCCProvider("RMSJouleHeatST");

//////////////////////////////////////////////////////////////////////////////

void RMSJouleHeatSourceTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("AddToBothEnergyEquations",
    "Add the RMS Joule heat term to both energy equations.");

  // Zuheyr for Argon radiative heat transfer source term
  options.addConfigOption< bool >(  "ComputeArgonRadiativeHeatTransferLoss","Argon Radiative heat transfer loss?");  
  options.addConfigOption< CFreal >("ArRadC1","Argon Radiative model polynomial constant.");
  options.addConfigOption< CFreal >("ArRadC2","Argon Radiative model polynomial constant.");
  options.addConfigOption< CFreal >("ArRadC3","Argon Radiative model polynomial constant.");
  options.addConfigOption< CFreal >("ArRadC4","Argon Radiative model polynomial constant.");
  options.addConfigOption< CFreal >("ArRadC5","Argon Radiative model polynomial constant.");
  options.addConfigOption< CFreal >("ArRadC6","Argon Radiative model polynomial constant.");
}

//////////////////////////////////////////////////////////////////////////////

RMSJouleHeatSourceTerm::RMSJouleHeatSourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  socket_rmsJouleHeatSource("rmsJouleHeatSource"),
  m_library(CFNULL),
  m_physicalData(),
  m_addToBothEnergyEquations(false),
  m_computeArgonRadiativeHeatTransferLoss(false),
  m_ArRadC1(0.0),
  m_ArRadC2(0.0),
  m_ArRadC3(0.0),
  m_ArRadC4(0.0),
  m_ArRadC5(0.0),
  m_ArRadC6(0.0)
{
 addConfigOptionsTo(this);

 m_addToBothEnergyEquations = false;
 setParameter("AddToBothEnergyEquations",&m_addToBothEnergyEquations);

 m_computeArgonRadiativeHeatTransferLoss = false;
 setParameter("ComputeArgonRadiativeHeatTransferLoss",&m_computeArgonRadiativeHeatTransferLoss);
 
 m_ArRadC1 = -6.932e-41;
 setParameter("ArRadC1",&m_ArRadC1);
 m_ArRadC2 =  4.753e-39;
 setParameter("ArRadC2",&m_ArRadC2);
 m_ArRadC3 = -5.808e-38;
 setParameter("ArRadC3",&m_ArRadC3);
 m_ArRadC4 = -2.843e-45;
 setParameter("ArRadC4",&m_ArRadC4);
 m_ArRadC5 =  8.629e-47;
 setParameter("ArRadC5",&m_ArRadC5);
 m_ArRadC6 =  6.706e+04;
 setParameter("ArRadC6",&m_ArRadC6);
}

//////////////////////////////////////////////////////////////////////////////

RMSJouleHeatSourceTerm::~RMSJouleHeatSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
RMSJouleHeatSourceTerm::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = ComputeSourceTermFVMCC::needsSockets();

  result.push_back(&socket_rmsJouleHeatSource);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void RMSJouleHeatSourceTerm::setup()
{
  using namespace COOLFluiD::Framework;

  ComputeSourceTermFVMCC::setup();
  
  m_library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(m_library.isNotNull());

  PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm()->resizePhysicalData(m_physicalData);
}

//////////////////////////////////////////////////////////////////////////////

void RMSJouleHeatSourceTerm::computeSource(Framework::GeometricEntity *const element,
					   RealVector& source,
					   RealMatrix& jacobian)
{
  using namespace COOLFluiD::Framework;
  
  // each cell, each iteration
  const CFreal qRad = computeRadiationLoss(element);  
  
  const CFuint nbEqs =
    PhysicalModelStack::getActive()->getEquationSubSysDescriptor().getNbEqsSS();
    
  // this is needed for coupling
  const CFuint nbEqsWoE = source.size() - 2;
  if (nbEqs == nbEqsWoE || nbEqs == source.size()) {
    CFLogDebugMin( "RMSJouleHeatSourceTerm::computeSource()" << "\n");
    
    const CFuint elemID = element->getID();
    DataHandle<CFreal> rmsJouleHeatSource = socket_rmsJouleHeatSource.getDataHandle();
    DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
    const CFuint nbTs = m_library->getNbTempVib() + m_library->getNbTe();
    const CFuint TID = source.size()-nbTs-2-1;
    // is not perturbed because it is computed in command, here is got just data handle
    if (nbTs == 0) {
     source[TID] = rmsJouleHeatSource[elemID] + qRad;
    }
    else {
      if (m_addToBothEnergyEquations) { 
        source[TID] = rmsJouleHeatSource[elemID] + qRad; 
      } 
      if (m_library->getNbTempVib() > 0 && m_library->getNbTe() == 0) {
	// we assume that Te == first Tv
	source[TID+1] = rmsJouleHeatSource[elemID] + qRad;
      }
      if (m_library->getNbTe() == 1) {
        source[TID+nbTs] = rmsJouleHeatSource[elemID] + qRad;
      }
    }
   
    CFLogDebugMax("RMSJouleHeatSourceTerm::computeSource() => source[" << TID << "] = " << source[TID] << "\n");
    
    source *= volumes[elemID];
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal RMSJouleHeatSourceTerm::computeRadiationLoss(Framework::GeometricEntity *const element)
{
  CFreal ArRadLoss = 0.;

  if (m_computeArgonRadiativeHeatTransferLoss) {

/*
 * CFout   << "Use: Simulator.SubSystem.CellCenterFVM.Data.RMSJouleHeatST.ComputeArgonRadiativeHeatTransferLoss = true/false \n";
 */  
    const CFuint nbSpecies = m_library->getNbSpecies();
    RealVector molarMasses(nbSpecies);
    RealVector ys(nbSpecies);
  
    m_library->getMolarMasses(molarMasses);
    m_library->getSpeciesMassFractions(ys);
  
    CFreal elY = ys[0];
    CFreal elM = molarMasses[0];
  
    Common::SafePtr<ConvectiveVarSet> updateVar = getMethodData().getUpdateVar();
    cf_assert(updateVar.isNotNull());
    Common::SafePtr<EulerTerm> eulerTerm = PhysicalModelStack::getActive()->getImplementor()->getConvectiveTerm().d_castTo<EulerTerm>();
    cf_assert(eulerTerm.isNotNull());
    const RealVector& refData = eulerTerm->getReferencePhysicalData();
    
    // set temperature and pressure
    State *currState =  element->getState(0);
    updateVar->computePhysicalData(*currState, m_physicalData);
  
    CFreal Tdim = m_physicalData[EulerTerm::T]*refData[EulerTerm::T];
    cf_assert(Tdim > 0.);
    cf_assert(Tdim < 15000.);
    CFreal pdim = eulerTerm->getPressureFromState(m_physicalData[EulerTerm::P])*refData[EulerTerm::P];
    cf_assert(pdim > 0.);
    CFreal ddim = m_physicalData[EulerTerm::RHO]*refData[EulerTerm::RHO];
    /* or   RealVector dhe(3);
     *      m_library->setDensityEnthalpyEnergy(Tdim,pdim,dhe);
     *      CFdouble density    = dhe[0];
     *      CFdouble enthalpyTt = dhe[1];
     */
    cf_assert(ddim > 0.);
    
    const CFdouble N_AVOGADRO = 6.022140857e23;
    CFreal elNumDen = N_AVOGADRO*elY*ddim/elM;
  
//  CFout << " option = m_computeArgonRadiativeHeatTransferLoss = " << m_computeArgonRadiativeHeatTransferLoss << "\n";
  /*
	Equation is:
	 ---------------------------------------------------------------------------
	| P=-(n_e)^2 (C1 + C2 * T^-0.5 + C3 * T^-1 + C4 * T + C5 * T^1.5)*exp(C6/T) |
	 ---------------------------------------------------------------------------
			ONLY FOR	3000 K < T < 15000 K
  */

    if (Tdim > 3000.) {
      ArRadLoss = - (elNumDen*elNumDen)*(m_ArRadC1 + m_ArRadC2*pow(Tdim,-0.5) + m_ArRadC3/Tdim + m_ArRadC4*Tdim + m_ArRadC5*pow(Tdim,1.5))*exp(m_ArRadC6/Tdim);
    }
  }

  return ArRadLoss;
}

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
