#include "FiniteVolumeNEQ/FiniteVolumeNEQ.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "FiniteVolumeNEQ/SubInletEuler2DProfileUVTYi.hh"
#include "Framework/MethodCommandProvider.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/MultiScalarTerm.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "NavierStokes/MultiScalarVarSet.hh"

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

MethodCommandProvider<SubInletEuler2DProfileUVTYi,
                      CellCenterFVMData,
                       FiniteVolumeNEQModule>
subInletEuler2DProfileUVTYiFVMCCProvider("SubInletEuler2DProfileUVTYi");

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler2DProfileUVTYi::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
  options.addConfigOption< std::vector<CFreal> >("Yi","mass fraction");
}

//////////////////////////////////////////////////////////////////////////////

SubInletEuler2DProfileUVTYi::SubInletEuler2DProfileUVTYi(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _library(CFNULL),
  _dataInnerState(),
  _bCoord(),
  _bState(3)
{
   addConfigOptionsTo(this);
  std::cout << "SubInletEuler2DProfileUVTYi" << std::endl;

  _functions = std::vector<std::string>();
  setParameter("Def",&_functions);

  _vars = std::vector<std::string>();
  setParameter("Vars",&_vars);

  _Yi = std::vector<CFreal>();
  setParameter("Yi",&_Yi);
}

//////////////////////////////////////////////////////////////////////////////

SubInletEuler2DProfileUVTYi::~SubInletEuler2DProfileUVTYi()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler2DProfileUVTYi::setGhostState(GeometricEntity *const face)
{
  State& innerState = *face->getState(0);
  State& ghostState = *face->getState(1);

  // coordinates of the boundary point
  _bCoord = (innerState.getCoordinates() + ghostState.getCoordinates());
  _bCoord *= 0.5;
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  cf_assert(nbDim == 2);
  for(CFuint iDim = 0; iDim < nbDim; iDim++){
    _variables[iDim] = _bCoord[iDim];
  }
  _variables[PhysicalModelStack::getActive()->getDim()] = SubSystemStatusStack::getActive()->getCurrentTimeDim();
  _vFunction.evaluate(_variables, _bState);

  _uinf = _bState[0];
  _vinf = _bState[1];
  _temperature = _bState[2];

  _uinf /= _varSet->getModel()->getVelRef();
  _vinf /= _varSet->getModel()->getVelRef();
  _temperature /= _varSet->getModel()->getTempRef();

  // set the physical data starting from the inner state
  _varSet->computePhysicalData(innerState, _dataInnerState);

  SafePtr<PhysicalChemicalLibrary> library = PhysicalModelStack::getActive()->getImplementor()->getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  assert(library.isNotNull());
  
  // physical constants
  const CFuint nbSpecies = _varSet->getModel()->getNbScalarVars(0);
  
  // impose velocity profile and temperature
  const CFreal uGhost = 2.0*_uinf - _dataInnerState[EulerTerm::VX];
  const CFreal vGhost = 2.0*_vinf - _dataInnerState[EulerTerm::VY];
  const CFreal TGhost = _temperature;
  
  // extrapolate pressure
  const CFreal pGhost = _dataInnerState[EulerTerm::P];
  
  library->getMolarMasses(_Mm);
  
  CFreal sumYiG = 0.;
  const CFuint firstSpecies = _varSet->getModel()->getFirstScalarVar(0);
  for (CFuint i = 0; i < nbSpecies; ++i){
    const CFreal yIGhost = 2.0*_Yi[i] - _dataInnerState[firstSpecies + i];
    sumYiG += yIGhost / _Mm[i];
    ghostState[i] = yIGhost;
  }   
  
  const CFreal rhoGhost =  pGhost/(library->getRgas()*TGhost*sumYiG);
  for (CFuint i = 0; i < nbSpecies; ++i){
    ghostState[i] *= rhoGhost; 
  }
  
  ghostState[nbSpecies]   = uGhost;
  ghostState[nbSpecies+1] = vGhost;
  ghostState[nbSpecies+2] = TGhost;
}

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler2DProfileUVTYi::setup()
{
  FVMCC_BC::setup();
   
  _bCoord.resize(PhysicalModelStack::getActive()->getDim());
  _variables.resize(PhysicalModelStack::getActive()->getDim() + 1);

  _varSet = getMethodData().getUpdateVar().d_castTo<MultiScalarVarSet<Euler2DVarSet> >();
  cf_assert(_varSet.isNotNull());
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
    
  //resize vector for U, V, T
  _bState.resize(3);
  _Yi.resize(_varSet->getModel()->getNbScalarVars(0));
  _Mm.resize(_varSet->getModel()->getNbScalarVars(0));
}

//////////////////////////////////////////////////////////////////////////////

void SubInletEuler2DProfileUVTYi::configure ( Config::ConfigArgs& args )
{
  FVMCC_BC::configure(args);

  _vFunction.setFunctions(_functions);
  _vFunction.setVariables(_vars);
  try {
    _vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
