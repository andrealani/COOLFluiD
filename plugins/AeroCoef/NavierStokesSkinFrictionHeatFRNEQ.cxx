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
}

//////////////////////////////////////////////////////////////////////////////

NavierStokesSkinFrictionHeatFRNEQ::NavierStokesSkinFrictionHeatFRNEQ(const std::string& name) :
  NavierStokesSkinFrictionHeatFluxFR(name),
  _tempVib()
{
  addConfigOptionsTo(this);
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
  Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  if (library->getNbTe() + library->getNbTempVib() > 0) {
    _tempVib.resize(library->getNbTe() + library->getNbTempVib());
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSkinFrictionHeatFRNEQ::computeDimensionalPressDensTemp
(CFreal& pDim, CFreal& rhoDim, CFreal& TDim, CFuint flxIdx)
{ 
  Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  const CFreal rhoRef = (m_updateVarSet->getModel()->getReferencePhysicalData())[EulerTerm::RHO];
  const CFuint nbSpecies = library->getNbSpecies();
  CFreal rho = 0.;
  for (CFuint i = 0; i < nbSpecies; ++i) {
    rho +=  (*m_cellStatesFlxPnt[flxIdx])[i];
  }
  rho *= rhoRef;
  
  TDim =  (*m_cellStatesFlxPnt[flxIdx])[m_TID] * (m_updateVarSet->getModel()->getTempRef());
  
  const CFuint startTID = this->m_TID + 1;
  for (CFuint i = 0; i <  _tempVib.size(); ++i) {
    _tempVib[i] =  (*m_cellStatesFlxPnt[flxIdx])[startTID + i]*(m_updateVarSet->getModel()->getTempRef());
  }
  
  CFreal* tVec = (_tempVib.size() == 0) ? CFNULL : &_tempVib[0];
  pDim = library->pressure(rho, TDim, tVec);
  rhoDim = m_rhoWall * rhoRef;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////




