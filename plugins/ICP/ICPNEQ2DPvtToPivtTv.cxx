#include "ICP/ICPNEQ.hh"
#include "ICPNEQ2DPvtToPivtTv.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ICP {

//////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ICPNEQ2DPvtToPivtTv, VarSetTransformer, ICPNEQModule,1>
icpNEQ2DPvtToPivtTvProvider("ICPNEQ2DPvtToPivtTv");

// AL: the following provider is defined for backward compatibility 
Environment::ObjectProvider<ICPNEQ2DPvtToPivtTv, VarSetTransformer, ICPNEQModule,1>
icpNEQ2DPuvtToPivtTvProvider("ICPNEQ2DPuvtToPivtTv");

//////////////////////////////////////////////////////////////////////

ICPNEQ2DPvtToPivtTv::ICPNEQ2DPvtToPivtTv
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  NEQ::Euler2DNEQPuvtToPivtTv(model)
{
}
      
//////////////////////////////////////////////////////////////////////
      
ICPNEQ2DPvtToPivtTv::~ICPNEQ2DPvtToPivtTv()
{
}

//////////////////////////////////////////////////////////////////////

void ICPNEQ2DPvtToPivtTv::transform(const State& state, State& result)
{
  NEQ::Euler2DNEQPuvtToPivtTv::transform(state, result);
  cf_assert(state.size() == result.size());
  
  const CFuint startE = PhysicalModelStack::getActive()->getDim()+2; 
  result[result.size()-2] = state[startE];
  result[result.size()-1] = state[startE+1];
}
      
//////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
