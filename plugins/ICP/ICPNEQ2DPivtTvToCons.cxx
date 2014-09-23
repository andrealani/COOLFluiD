#include "ICP/ICPNEQ.hh"
#include "ICPNEQ2DPivtTvToCons.hh"
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

    namespace ICP {

//////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ICPNEQ2DPivtTvToCons, VarSetTransformer, ICPNEQModule,1>
icpNEQ2DPivtTvToConsProvider("ICPNEQ2DPivtTvToCons");

//////////////////////////////////////////////////////////////////////

ICPNEQ2DPivtTvToCons::ICPNEQ2DPivtTvToCons
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  Euler2DNEQPivtTvToCons(model)
{
  cf_assert(model.isNotNull());
}
      
//////////////////////////////////////////////////////////////////////
      
ICPNEQ2DPivtTvToCons::~ICPNEQ2DPivtTvToCons()
{
}

//////////////////////////////////////////////////////////////////////

void ICPNEQ2DPivtTvToCons::transform(const State& state, State& result)
{
  Euler2DNEQPivtTvToCons::transform(state, result);

  // here we assume that E's are the last two components
  const CFuint startE = result.size()-2;
  result[startE]   = state[startE];
  result[startE+1] = state[startE+1];
}
      
//////////////////////////////////////////////////////////////////////
      
void ICPNEQ2DPivtTvToCons::transformFromRef(const RealVector& data, State& result)
{ 
  Euler2DNEQPivtTvToCons::transformFromRef(data, result);

  // here we assume that E's are the last two components
  const CFuint firstScalarVar = _model->getDataSize() - 2;
  const CFuint startE = result.size()-2;
  cf_assert(firstScalarVar < data.size());
  result[startE]   = data[firstScalarVar];
  result[startE+1] = data[firstScalarVar+1];
}

//////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
