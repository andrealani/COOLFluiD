#include "Common/PE.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/DirPaths.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/BadValueException.hh"
#include "Framework/PathAppender.hh"
#include "FiniteVolume/CellCenterFVM.hh"
#include "FiniteVolume/DerivativeComputer.hh"
#include "NavierStokes/NavierStokesVarSet.hh"
#include "NEQ/NavierStokesNEQVarSet.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "AeroCoef/AeroCoefFVM.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "AeroCoef/NavierStokesSkinFrictionHeatFluxCCNEQ.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Numerics::FiniteVolume;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NavierStokesSkinFrictionHeatFluxCCNEQ,
          DataProcessingData,
          AeroCoefFVMModule>
navierStokesSkinFrictionHeatFluxCCNEQProvider
("NavierStokesSkinFrictionHeatFluxCCNEQ");

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSkinFrictionHeatFluxCCNEQ::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

NavierStokesSkinFrictionHeatFluxCCNEQ::NavierStokesSkinFrictionHeatFluxCCNEQ(const std::string& name) :
  NavierStokesSkinFrictionHeatFluxCC(name),
  _tempVib()
{
  addConfigOptionsTo(this);
}
      
//////////////////////////////////////////////////////////////////////////////

NavierStokesSkinFrictionHeatFluxCCNEQ::~NavierStokesSkinFrictionHeatFluxCCNEQ()
{
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSkinFrictionHeatFluxCCNEQ::setup()
{
  CFAUTOTRACE;

  NavierStokesSkinFrictionHeatFluxCC::setup();
  Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  if (library->getNbTe() + library->getNbTempVib() > 0) {
    _tempVib.resize(library->getNbTe() + library->getNbTempVib());
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesSkinFrictionHeatFluxCCNEQ::computeDimensionalPressDensTemp
(CFreal& pDim, CFreal& rhoDim, CFreal& TDim)
{ 
  Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  const CFreal rhoRef = (m_updateVarSet->getModel()->getReferencePhysicalData())[EulerTerm::RHO];
  const CFuint nbSpecies = library->getNbSpecies();
  CFreal rho = 0.;
  for (CFuint i = 0; i < nbSpecies; ++i) {
    rho += (*m_avState)[i];
  }
  rho *= rhoRef;
  
  TDim = (*m_avState)[m_TID] * (m_updateVarSet->getModel()->getTempRef());
  
  const CFuint startTID = this->m_TID + 1;
  for (CFuint i = 0; i <  _tempVib.size(); ++i) {
    _tempVib[i] = (*m_avState)[startTID + i]*(m_updateVarSet->getModel()->getTempRef());
  }
  
  CFreal* tVec = (_tempVib.size() == 0) ? CFNULL : &_tempVib[0];
  pDim = library->pressure(rho, TDim, tVec);
  rhoDim = _rhoWall * rhoRef;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////




