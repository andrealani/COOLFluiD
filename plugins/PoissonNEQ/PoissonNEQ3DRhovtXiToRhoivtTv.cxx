#include "PoissonNEQ/PoissonNEQ.hh"
#include "PoissonNEQ3DRhovtXiToRhoivtTv.hh"
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

Environment::ObjectProvider<PoissonNEQ3DRhovtXiToRhoivtTv, VarSetTransformer, PoissonNEQModule,1>
PoissonNEQ3DRhovtXiToRhoivtTvProvider("PoissonNEQ3DRhovtXiToRhoivtTv");

//////////////////////////////////////////////////////////////////////

PoissonNEQ3DRhovtXiToRhoivtTv::PoissonNEQ3DRhovtXiToRhoivtTv
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  Euler3DNEQRhovtXiToRhoivtTv(model)
{
  cf_assert(model.isNotNull());
}
      
//////////////////////////////////////////////////////////////////////
      
PoissonNEQ3DRhovtXiToRhoivtTv::~PoissonNEQ3DRhovtXiToRhoivtTv()
{
}

//////////////////////////////////////////////////////////////////////

void PoissonNEQ3DRhovtXiToRhoivtTv::transform(const State& state, State& result)
{
  //cout<<"PoissonNEQ3DRhovtXiToRhoivtTv\n";
  Euler3DNEQRhovtXiToRhoivtTv::transform(state, result);
  //cout<<"PoissonNEQ3DRhovtXiToRhoivtTv2\n";
  // here we assume that phi is the last component
  const CFuint startE = result.size()-1;
  //CFLog(INFO, "startE = " <<startE <<"\n");
  //CFLog(INFO, "state[startE] = " <<state[startE] <<"\n");
  result[startE]   = state[startE];
}
      
//////////////////////////////////////////////////////////////////////
      
void PoissonNEQ3DRhovtXiToRhoivtTv::transformFromRef(const RealVector& data, State& result)
{ 
  //cout<<"here 8b\n";
  Euler3DNEQRhovtXiToRhoivtTv::transformFromRef(data, result);

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
