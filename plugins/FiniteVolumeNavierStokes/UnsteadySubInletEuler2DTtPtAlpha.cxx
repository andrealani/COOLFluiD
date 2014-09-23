#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "FiniteVolumeNavierStokes/UnsteadySubInletEuler2DTtPtAlpha.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"

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

MethodCommandProvider<UnsteadySubInletEuler2DTtPtAlpha, CellCenterFVMData, 
FiniteVolumeNavierStokesModule> UnsteadySubInletEuler2DTtPtAlphaFVMCCProvider("UnsteadySubInletEuler2DTtPtAlphaFVMCC");

//////////////////////////////////////////////////////////////////////////////

void UnsteadySubInletEuler2DTtPtAlpha::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
}

//////////////////////////////////////////////////////////////////////////////

UnsteadySubInletEuler2DTtPtAlpha::UnsteadySubInletEuler2DTtPtAlpha(const std::string& name) :
  FVMCC_BC(name),
  _bCoord(),
  _variables(),
  _values(),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
   addConfigOptionsTo(this);

  _functions = std::vector<std::string>();
   setParameter("Def",&_functions);

  _vars = std::vector<std::string>();
   setParameter("Vars",&_vars);
}

//////////////////////////////////////////////////////////////////////////////

UnsteadySubInletEuler2DTtPtAlpha::~UnsteadySubInletEuler2DTtPtAlpha()
{
}

//////////////////////////////////////////////////////////////////////////////

void UnsteadySubInletEuler2DTtPtAlpha::configure ( Config::ConfigArgs& args )
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

void UnsteadySubInletEuler2DTtPtAlpha::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  for(CFuint iDim = 0; iDim < nbDim; iDim++)
  {
    _variables[iDim] = _bCoord[iDim];
  }
  _variables[PhysicalModelStack::getActive()->getDim()] = SubSystemStatusStack::getActive()->getCurrentTimeDim();

  //Evaluate the function
  _vFunction.evaluate(_variables,_values);

  _tTotal = _values[0];
  _pTotal = _values[1];
  _alpha = _values[2];

  _pTotal/=_varSet->getModel()->getPressRef();
  _tTotal/=_varSet->getModel()->getTempRef();


  // set the physical data starting from the inner state
  _varSet->computePhysicalData(*innerState, _dataInnerState);

  const CFreal gamma = _varSet->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal gammaDivGammaMinus1 = gamma/gammaMinus1;
  const CFreal R = _varSet->getModel()->getR();
  const CFreal u = _dataInnerState[EulerTerm::VX];
  const CFreal v = _dataInnerState[EulerTerm::VY];
  const CFreal uSqvSq = u*u + v*v;
  const CFreal pInnerState = _dataInnerState[EulerTerm::P];
  const CFreal tInnerState = pInnerState / (R*_dataInnerState[EulerTerm::RHO]);
  const CFreal machInner = sqrt(uSqvSq/(gamma*R*tInnerState));
  const CFreal coefficient = 1.0 + 0.5*gammaMinus1*machInner*machInner;
  const CFreal tTotalInner = tInnerState*coefficient;
  const CFreal coeffPow = pow(coefficient, gamma/gammaMinus1);
  const CFreal pTotalInner = pInnerState*coeffPow;
  const CFreal tgAlphaInner = v/u;

  // ghost state quantities
  const CFreal tTotalGhost = 2.0*_tTotal - tTotalInner;
  const CFreal pTotalGhost = 2.0*_pTotal - pTotalInner;
  const CFreal tgAlphaGhost = 2.0*tan(_alpha) - tgAlphaInner;
  const CFreal machGhost = machInner;
  const CFreal tGhost = tTotalGhost / coefficient;

  // set the physical data for the ghost state
  _dataGhostState[EulerTerm::P] = pTotalGhost/coeffPow;
  _dataGhostState[EulerTerm::RHO] = _dataGhostState[EulerTerm::P]/(R*tGhost);
  _dataGhostState[EulerTerm::VX] = machGhost*sqrt(gamma*R*tGhost/
						  (1.0 + tgAlphaGhost*tgAlphaGhost));
  _dataGhostState[EulerTerm::VY] = tgAlphaGhost*_dataGhostState[EulerTerm::VX];
  _dataGhostState[EulerTerm::V] = sqrt(_dataGhostState[EulerTerm::VX]*
				 _dataGhostState[EulerTerm::VX] +
				 _dataGhostState[EulerTerm::VY]*
				 _dataGhostState[EulerTerm::VY]);
  _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*_dataGhostState[EulerTerm::P]
				   + 0.5*_dataGhostState[EulerTerm::RHO]*
				   _dataGhostState[EulerTerm::V]*
				   _dataGhostState[EulerTerm::V])/
    _dataGhostState[EulerTerm::RHO];
  _dataGhostState[EulerTerm::A] = sqrt(gamma*_dataGhostState[EulerTerm::P]/
				       _dataGhostState[EulerTerm::RHO]);
  
  _dataGhostState[EulerTerm::T] = tGhost;
  
  // set the ghost state starting from the physical data
  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
 }

//////////////////////////////////////////////////////////////////////////////

void UnsteadySubInletEuler2DTtPtAlpha::setup()
{
  FVMCC_BC::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);

  _bCoord.resize(PhysicalModelStack::getActive()->getDim());
  _variables.resize(PhysicalModelStack::getActive()->getDim() + 1);

  //resize vector for Tt, pt and alpha
  _values.resize(3);

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
