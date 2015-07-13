#include "NavierStokes/NavierStokes.hh"
#include "NavierStokes/EulerPvtToCons.hh"
#include "NavierStokes/EulerTerm.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<EulerPvtToCons, VarSetTransformer, NavierStokesModule, 1> 
euler1DPvtToConsProvider("Euler1DPvtToCons");

Environment::ObjectProvider<EulerPvtToCons, VarSetTransformer, NavierStokesModule, 1> 
euler2DPuvtToConsProvider("Euler2DPuvtToCons");

Environment::ObjectProvider<EulerPvtToCons, VarSetTransformer, NavierStokesModule, 1> 
euler3DPvtToConsProvider("Euler3DPvtToCons");

//////////////////////////////////////////////////////////////////////////////

EulerPvtToCons::EulerPvtToCons(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>()),
  m_rho(0.)
{
}

//////////////////////////////////////////////////////////////////////////////

EulerPvtToCons::~EulerPvtToCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void EulerPvtToCons::transform(const State& state, State& result)
{
  const CFreal p = _model->getPressureFromState(state[0]);
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint TID = dim+1;
  const CFreal T = state[TID];
  m_rho = _model->getDensity(p,T);
  
  result[0] = m_rho;
  
  CFreal V2 = 0.;
  for (CFuint i=0; i < dim; ++i) {
    const CFuint uID = i+1;
    const CFreal u = state[uID];
    result[uID] = m_rho*u; 
    V2 += u*u;
  }
  
  result[TID] = p/(_model->getGamma() - 1.) + 0.5*m_rho*V2;
}
      
//////////////////////////////////////////////////////////////////////////////

void EulerPvtToCons::transformFromRef(const RealVector& data, State& result)
{
  m_rho = data[EulerTerm::RHO];

  result[0] = m_rho;
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  for (CFuint i=0; i < dim; ++i) {
    result[i+1] = m_rho*data[EulerTerm::VX+i];
  }
  
  const CFuint TID = dim+1;
  result[TID] = m_rho*data[EulerTerm::H] - _model->getPressureFromState(data[EulerTerm::P]);
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
