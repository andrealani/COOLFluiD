#include "NEQ/NEQ.hh"
#include "Euler2DNEQMachAlphaPTyiToRhoivtTv.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DNEQMachAlphaPTyiToRhoivtTv, VarSetTransformer, NEQModule, 1> 
euler2DNEQMachAlphaPTyiToRhoivtTvProvider("Euler2DNEQMachAlphaPTyiToRhoivtTv");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQMachAlphaPTyiToRhoivtTv::Euler2DNEQMachAlphaPTyiToRhoivtTv
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<NEQTerm>()),
  _tvDim(),
  _mmasses(),
  _ys()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQMachAlphaPTyiToRhoivtTv::~Euler2DNEQMachAlphaPTyiToRhoivtTv()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQMachAlphaPTyiToRhoivtTv::setup(const CFuint maxNbTransStates)
{
  VarSetTransformer::setup(maxNbTransStates);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQMachAlphaPTyiToRhoivtTv::transform(const State& state, State& result)
{ 
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();

  const CFuint nbSpecies = _model->getNbScalarVars(0);
  
  /// AL: wild approximation here: we consider the isoentropic relations valid also for
  /// thermo-chemical non equilibrium!!!
  const CFreal Mach = state[0];
  // convert angle in radiants
  const CFreal alpha = state[1]*MathTools::MathConsts::CFrealPi()/180.;
  
  _mmasses.resize(nbSpecies);
  _ys.resize(nbSpecies);
  
  library->getMolarMasses(_mmasses);
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbY = nbSpecies - 1;
  const CFuint startY = nbEqs - nbY;
  
  // Set the mixture density (sum of the partial densities)
  CFreal invM = 0.0;
  CFreal sumy = 0.0;
  for (CFuint is = 0; is < nbY; ++is) {
    cf_assert(state[startY + is] < 1.00001);
    invM += state[startY + is]/_mmasses[is];
    sumy += state[startY + is];
    _ys[is] = state[startY + is];
  }
  // take into account the last mass fraction
  invM += (1. - sumy)/_mmasses[nbY];
  _ys[nbY]  = (1. - sumy);
  
  const CFreal Rgas = library->getRgas();
  CFreal p = state[2];
  CFreal T = state[3]; 
  CFreal rho = p/(Rgas*T*invM);
  
  for (CFuint is = 0; is < nbY; ++is) {
    result[is] = rho*state[startY + is];
  }
  // take into account the last mass fraction
  result[nbY] = rho*(1. - sumy);
  
  CFreal gamma = 0.0;
  CFreal a = 0.0;
  
  // set the vibrational temperature
  const CFuint nbTv = _model->getNbScalarVars(1);
  _tvDim.resize(nbTv);
  for (CFuint i = 0; i < nbTv; ++i) {
    result[nbSpecies + 3 + i] = _tvDim[i] = state[4 + i];
  }
  
  // set the current species fractions in the thermodynamic library
  // this has to be done right here, before computing any other thermodynamic quantity !!! 
  library->setSpeciesFractions(_ys);
  
  // @TODO/ here adapt for adimensionalization 
  library->frozenGammaAndSoundSpeed(T,p,rho,gamma, a, &_tvDim);
  const CFreal u = cos(alpha)*(Mach * a);
  const CFreal v = sin(alpha)*(Mach * a);
    
  result[nbSpecies] = u;
  result[nbSpecies + 1] = v;
  result[nbSpecies + 2] = T;
}
      
//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQMachAlphaPTyiToRhoivtTv::transformFromRef(const RealVector& data,State& result)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DNEQMachAlphaPTyiToRhoivtTv::transformFromRef()");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
