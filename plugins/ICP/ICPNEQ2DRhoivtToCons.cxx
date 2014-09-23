#include "ICP/ICPNEQ.hh"
#include "ICPNEQ2DRhoivtToCons.hh"
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

Environment::ObjectProvider<ICPNEQ2DRhoivtToCons, VarSetTransformer, ICPNEQModule,1>
icpNEQ2DRhoivtToConsProvider("ICPNEQ2DRhoivtToCons");

//////////////////////////////////////////////////////////////////////

ICPNEQ2DRhoivtToCons::ICPNEQ2DRhoivtToCons
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  Euler2DNEQRhoivtToCons(model)
{
  cf_assert(model.isNotNull());
}
      
//////////////////////////////////////////////////////////////////////
      
ICPNEQ2DRhoivtToCons::~ICPNEQ2DRhoivtToCons()
{
}

//////////////////////////////////////////////////////////////////////

void ICPNEQ2DRhoivtToCons::transform(const State& state, State& result)
{
  Euler2DNEQRhoivtToCons::transform(state, result);

  // here we assume that E's are the last two components
  const CFuint startE = result.size()-2;
  result[startE]   = state[startE];
  result[startE+1] = state[startE+1];
}
      
//////////////////////////////////////////////////////////////////////
      
void ICPNEQ2DRhoivtToCons::transformFromRef(const RealVector& data, State& result)
{ 
  Euler2DNEQRhoivtToCons::transformFromRef(data, result);

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
