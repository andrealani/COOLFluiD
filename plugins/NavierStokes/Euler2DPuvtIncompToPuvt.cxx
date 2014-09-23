#include "NavierStokes/NavierStokes.hh"
#include "Euler2DPuvtIncompToPuvt.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NavierStokes {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DPuvtIncompToPuvt, VarSetTransformer, NavierStokesModule, 1> 
euler2DPuvtIncompToPuvtProvider("Euler2DPuvtIncompToPuvt");

//////////////////////////////////////////////////////////////////////////////

Euler2DPuvtIncompToPuvt::Euler2DPuvtIncompToPuvt(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DPuvtIncompToPuvt::~Euler2DPuvtIncompToPuvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvtIncompToPuvt::transform(const State& state, State& result)
{
  result[0] = state[0] + 10000.; // @TODO AL: find a better way!!!
  result[1] = state[1];
  result[2] = state[2];
  result[3] = state[3];
  std::cout << "################################ ciao #########"<< std::endl; 
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DPuvtIncompToPuvt::transformFromRef(const RealVector& data, State& result)
{
  result[0] = data[EulerTerm::P] + 10000.; // @TODO AL: find a better way!!!
  result[1] = data[EulerTerm::VX];
  result[2] = data[EulerTerm::VY];
  result[3] = data[EulerTerm::T];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
