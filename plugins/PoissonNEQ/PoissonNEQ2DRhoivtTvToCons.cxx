#include "PoissonNEQ/PoissonNEQ.hh"
#include "PoissonNEQ2DRhoivtTvToCons.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::NEQ;

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace PoissonNEQ {

//////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<PoissonNEQ2DRhoivtTvToCons, VarSetTransformer, PoissonNEQModule,1>
poissonNEQ2DRhoivtTvToConsProvider("PoissonNEQ2DRhoivtTvToCons");

//////////////////////////////////////////////////////////////////////

PoissonNEQ2DRhoivtTvToCons::PoissonNEQ2DRhoivtTvToCons
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  Euler2DNEQRhoivtTvToCons(model)
{
  cf_assert(model.isNotNull());
}
      
//////////////////////////////////////////////////////////////////////
      
PoissonNEQ2DRhoivtTvToCons::~PoissonNEQ2DRhoivtTvToCons()
{
}

//////////////////////////////////////////////////////////////////////

void PoissonNEQ2DRhoivtTvToCons::transform(const State& state, State& result)
{
  Euler2DNEQRhoivtTvToCons::transform(state, result);

  // here we assume that phi is the last component
  const CFuint startE = result.size()-1;
  result[startE]   = state[startE];
}
      
//////////////////////////////////////////////////////////////////////
      
void PoissonNEQ2DRhoivtTvToCons::transformFromRef(const RealVector& data, State& result)
{ 
  Euler2DNEQRhoivtTvToCons::transformFromRef(data, result);

  // here we assume that phi is the last component
  const CFuint firstScalarVar = _model->getDataSize() - 1;
  const CFuint startE = result.size()-1;
  cf_assert(firstScalarVar < data.size());
  result[startE]   = data[firstScalarVar];
}

//////////////////////////////////////////////////////////////////////

    } // namespace PoissonNEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
