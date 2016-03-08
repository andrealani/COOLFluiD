#include "ICP/ICPNEQ.hh"
#include "ICP/ICPNEQ2DRhoivtTvToPivtTv.hh"
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

Environment::ObjectProvider<ICPNEQ2DRhoivtTvToPivtTv, VarSetTransformer, ICPNEQModule,1>
icpNEQ2DRhoivtTvToPivtTvProvider("ICPNEQ2DRhoivtTvToPivtTv");

//////////////////////////////////////////////////////////////////////

ICPNEQ2DRhoivtTvToPivtTv::ICPNEQ2DRhoivtTvToPivtTv
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  Euler2DNEQRhoivtTvToPivtTv(model)
{
  cf_assert(model.isNotNull());
}
      
//////////////////////////////////////////////////////////////////////
      
ICPNEQ2DRhoivtTvToPivtTv::~ICPNEQ2DRhoivtTvToPivtTv()
{
}

//////////////////////////////////////////////////////////////////////

void ICPNEQ2DRhoivtTvToPivtTv::transform(const State& state, State& result)
{
  Euler2DNEQRhoivtTvToPivtTv::transform(state, result);
  
  // here we assume that E's are the last two components
  const CFuint startE = result.size()-2;
  result[startE]   = state[startE];
  result[startE+1] = state[startE+1];
  
  CFLog(DEBUG_MAX, "ICPNEQ2DRhoivtTvToPivtTv::transform() => AFTER result = [" << result << "]\n");
}
      
//////////////////////////////////////////////////////////////////////
      
void ICPNEQ2DRhoivtTvToPivtTv::transformFromRef(const RealVector& data, State& result)
{ 
  Euler2DNEQRhoivtTvToPivtTv::transformFromRef(data, result);
  
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
