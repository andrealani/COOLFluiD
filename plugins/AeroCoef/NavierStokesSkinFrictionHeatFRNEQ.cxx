#include <algorithm>
#include <limits>
#include <cmath>

#include "Common/PE.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/DirPaths.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/BadValueException.hh"
#include "Framework/PathAppender.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolver.hh"
#include "NavierStokes/NavierStokesVarSet.hh"
#include "NEQ/NavierStokesNEQVarSet.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "AeroCoef/AeroCoefFR.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "AeroCoef/NavierStokesSkinFrictionHeatFRNEQ.hh"
#include "AeroCoef/NavierStokesSkinFrictionHeatFluxFR.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
//using namespace COOLFluiD::Numerics::FiniteVolume;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NavierStokesSkinFrictionHeatFRNEQ,
          DataProcessingData,
          AeroCoefFRModule>
navierStokesSkinFrictionHeatFRNEQProvider
("NavierStokesSkinFrictionHeatFRNEQ");

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSkinFrictionHeatFRNEQ::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption<CFreal>("TotalEnthalpyInf","Freestream total enthalpy [J/kg], including chemistry and modes, required for NEQ StantonNumberID=1.");
}

//////////////////////////////////////////////////////////////////////////////

NavierStokesSkinFrictionHeatFRNEQ::NavierStokesSkinFrictionHeatFRNEQ(const std::string& name) :
  NavierStokesSkinFrictionHeatFluxFR(name),
  _tempVib()
{
  addConfigOptionsTo(this);

  m_totalEnthalpyInf = std::numeric_limits< CFreal >::quiet_NaN();
  setParameter("TotalEnthalpyInf",&m_totalEnthalpyInf);
}
      
//////////////////////////////////////////////////////////////////////////////

NavierStokesSkinFrictionHeatFRNEQ::~NavierStokesSkinFrictionHeatFRNEQ()
{
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSkinFrictionHeatFRNEQ::setup()
{
  CFAUTOTRACE;

  NavierStokesSkinFrictionHeatFluxFR::setup();

  if (PhysicalModelStack::getActive()->getImplementor()->isAdimensional() ||
      PhysicalModelStack::getActive()->getImplementor()->getRefLength() != 1.)
  {
    throw BadValueException(FromHere(),"NavierStokesSkinFrictionHeatFRNEQ: adimensional variables or refLength != 1 are not supported; use dimensional variables with refLength = 1.");
  }

  if (m_stantonNumID != 0 && m_stantonNumID != 1)
  {
    throw BadValueException(FromHere(),"NavierStokesSkinFrictionHeatFRNEQ: this StantonNumberID is not supported; use 0 or 1.");
  }

  if (m_stantonNumID == 1 && !std::isfinite(m_totalEnthalpyInf))
  {
    throw BadValueException(FromHere(),"NavierStokesSkinFrictionHeatFRNEQ: StantonNumberID = 1 without TotalEnthalpyInf is not supported; set TotalEnthalpyInf.");
  }

  Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  if (library->getNbTe() + library->getNbTempVib() > 0) {
    _tempVib.resize(library->getNbTe() + library->getNbTempVib());
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal NavierStokesSkinFrictionHeatFRNEQ::computeStantonNumber(CFreal heatFlux, CFreal temperature, CFuint flxIdx)
{
  if (m_stantonNumID == 0)
  {
    return heatFlux/(m_rhoInf*std::pow(m_uInf,3.));
  }

  // static enthalpy at the wall, h_w = H - 0.5 |u|^2
  m_updateVarSet->computePhysicalData(*m_cellStatesFlxPnt[flxIdx],m_dataState);

  CFreal kineticEnergy = 0.;
  for (CFuint iDim = 0; iDim < m_dim; ++iDim)
  {
    kineticEnergy += 0.5*m_dataState[EulerTerm::VX+iDim]*m_dataState[EulerTerm::VX+iDim];
  }

  const CFreal wallEnthalpy = m_dataState[EulerTerm::H] - kineticEnergy;

  return heatFlux/(m_rhoInf*m_uInf*(m_totalEnthalpyInf-wallEnthalpy));
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSkinFrictionHeatFRNEQ::computeDimensionalPressDensTemp
(CFreal& pDim, CFreal& rhoDim, CFreal& TDim, CFuint flxIdx)
{
  m_diffVar->setComposition(*m_cellStatesFlxPnt[flxIdx],false,0);

  Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  const CFreal rhoRef = (m_updateVarSet->getModel()->getReferencePhysicalData())[EulerTerm::RHO];
  // the diffusive var set knows how the species are stored (rho_i or ln rho_i)
  CFreal rho = m_diffVar->getDensity(*m_cellStatesFlxPnt[flxIdx])*rhoRef;
  
  // logarithmic update variables (LogRhoivLogTTv: ln rho_i, u, v, ln T, ln Tv) store ln T and ln Tv
  const std::vector<std::string>& varNames = m_updateVarSet->getVarNames();
  const bool logVariables = (std::count(varNames.begin(), varNames.end(), "lnrho0") > 0);
  const RealVector& state = *m_cellStatesFlxPnt[flxIdx];
  const CFreal TRef = m_updateVarSet->getModel()->getTempRef();

  TDim = (logVariables ? std::exp(state[m_TID]) : state[m_TID])*TRef;
  
  const CFuint startTID = this->m_TID + 1;
  for (CFuint i = 0; i <  _tempVib.size(); ++i) {
    const CFreal Tv = state[startTID + i];
    _tempVib[i] = (logVariables ? std::exp(Tv) : Tv)*TRef;
  }
  
  CFreal* tVec = (_tempVib.size() == 0) ? CFNULL : &_tempVib[0];
  pDim = library->pressure(rho, TDim, tVec);
  rhoDim = m_rhoWall * rhoRef;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////




