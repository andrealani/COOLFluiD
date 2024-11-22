#include "PoissonNEQ/PoissonNEQ.hh"
#include "PoissonNEQ3DRhoivtToCons.hh"
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

Environment::ObjectProvider<PoissonNEQ3DRhoivtToCons, VarSetTransformer, PoissonNEQModule,1>
PoissonNEQ3DRhoivtToConsProvider("PoissonNEQ3DRhoivtToCons");

//////////////////////////////////////////////////////////////////////

PoissonNEQ3DRhoivtToCons::PoissonNEQ3DRhoivtToCons
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  Euler3DNEQRhoivtToCons(model)
{
  cf_assert(model.isNotNull());
}
      
//////////////////////////////////////////////////////////////////////
      
PoissonNEQ3DRhoivtToCons::~PoissonNEQ3DRhoivtToCons()
{
}

//////////////////////////////////////////////////////////////////////

void PoissonNEQ3DRhoivtToCons::transform(const State& state, State& result)
{
  //cout<< "PoissonNEQ3DRhoivtToCons__\n";
  Euler3DNEQRhoivtToCons::transform(state, result);
  //cout<< "PoissonNEQ3DRhoivtToCons\n";
  // here we assume that phi is the last component
  const CFuint startE = result.size()-1;
  result[startE]   = state[startE];
}
      
//////////////////////////////////////////////////////////////////////
      
void PoissonNEQ3DRhoivtToCons::transformFromRef(const RealVector& data, State& result)
{ 
  Euler3DNEQRhoivtToCons::transformFromRef(data, result);

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
