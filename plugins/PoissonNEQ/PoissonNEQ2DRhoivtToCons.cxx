#include "PoissonNEQ/PoissonNEQ.hh"
#include "PoissonNEQ2DRhoivtToCons.hh"
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

Environment::ObjectProvider<PoissonNEQ2DRhoivtToCons, VarSetTransformer, PoissonNEQModule,1>
PoissonNEQ2DRhoivtToConsProvider("PoissonNEQ2DRhoivtToCons");

//////////////////////////////////////////////////////////////////////

PoissonNEQ2DRhoivtToCons::PoissonNEQ2DRhoivtToCons
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  Euler2DNEQRhoivtToCons(model)
{
  cf_assert(model.isNotNull());
}
      
//////////////////////////////////////////////////////////////////////
      
PoissonNEQ2DRhoivtToCons::~PoissonNEQ2DRhoivtToCons()
{
}

//////////////////////////////////////////////////////////////////////

void PoissonNEQ2DRhoivtToCons::transform(const State& state, State& result)
{
  //cout<< "PoissonNEQ2DRhoivtToCons__\n";
  Euler2DNEQRhoivtToCons::transform(state, result);
 // cout<< "PoissonNEQ2DRhoivtToCons\n";
  // here we assume that phi is the last component
  const CFuint startE = result.size()-1;
  result[startE]   = state[startE];
}
      
//////////////////////////////////////////////////////////////////////
      
void PoissonNEQ2DRhoivtToCons::transformFromRef(const RealVector& data, State& result)
{ 
  Euler2DNEQRhoivtToCons::transformFromRef(data, result);

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
