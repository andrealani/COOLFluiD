#include "ICP/ICPNEQ.hh"
#include "ICPNEQ2DPuvtToRhoivtTv.hh"
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

Environment::ObjectProvider<ICPNEQ2DPuvtToRhoivtTv, VarSetTransformer, ICPNEQModule,1>
icpNEQ2DPuvtToRhoivtTvProvider("ICPNEQ2DPuvtToRhoivtTv");

//////////////////////////////////////////////////////////////////////

ICPNEQ2DPuvtToRhoivtTv::ICPNEQ2DPuvtToRhoivtTv
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  NEQ::Euler2DNEQPuvtToRhoivtTv(model)
{
}
      
//////////////////////////////////////////////////////////////////////
      
ICPNEQ2DPuvtToRhoivtTv::~ICPNEQ2DPuvtToRhoivtTv()
{
}

//////////////////////////////////////////////////////////////////////

void ICPNEQ2DPuvtToRhoivtTv::transform(const State& state, State& result)
{
  NEQ::Euler2DNEQPuvtToRhoivtTv::transform(state, result);
  cf_assert(state.size() == result.size());
  
  const CFuint startE = PhysicalModelStack::getActive()->getDim()+2; 
  result[result.size()-2] = state[startE];
  result[result.size()-1] = state[startE+1];
  
  // cout << "result ER EI = " <<  result[result.size()-2]  << " " << result[result.size()-1]  << endl;
  // cout << "state = " <<  state << endl;
}
      
//////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
