#include "PoissonNEQ/PoissonNEQ.hh"
#include "PoissonNEQ3DRhoivtTvToCons.hh"
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

Environment::ObjectProvider<PoissonNEQ3DRhoivtTvToCons, VarSetTransformer, PoissonNEQModule,1>
PoissonNEQ3DRhoivtTvToConsProvider("PoissonNEQ3DRhoivtTvToCons");

//////////////////////////////////////////////////////////////////////

PoissonNEQ3DRhoivtTvToCons::PoissonNEQ3DRhoivtTvToCons
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  Euler3DNEQRhoivtTvToCons(model)
{
  cf_assert(model.isNotNull());
}
      
//////////////////////////////////////////////////////////////////////
      
PoissonNEQ3DRhoivtTvToCons::~PoissonNEQ3DRhoivtTvToCons()
{
}

//////////////////////////////////////////////////////////////////////

void PoissonNEQ3DRhoivtTvToCons::transform(const State& state, State& result)
{
  //cout<<"PoissonNEQ3DRhoivtTvToCons\n";
  Euler3DNEQRhoivtTvToCons::transform(state, result);
  //cout<<"PoissonNEQ3DRhoivtTvToCons2\n";
  // here we assume that phi is the last component
  const CFuint startE = result.size()-1;
  result[startE]   = state[startE];
}
      
//////////////////////////////////////////////////////////////////////
      
void PoissonNEQ3DRhoivtTvToCons::transformFromRef(const RealVector& data, State& result)
{ 
  //cout<<"here 8b\n";
  Euler3DNEQRhoivtTvToCons::transformFromRef(data, result);

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
