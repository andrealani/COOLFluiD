#include "ICP/ICP.hh"
#include "ICPLTE2DPuvtToCons.hh"
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

Environment::ObjectProvider<ICPLTE2DPuvtToCons, VarSetTransformer, ICPModule,1>
icpLTE2DPuvtToConsProvider("ICPLTE2DPuvtToCons");

//////////////////////////////////////////////////////////////////////

ICPLTE2DPuvtToCons::ICPLTE2DPuvtToCons
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _icpModel(model->getConvectiveTerm().d_castTo<PTERM>()),
  _dhe(3)
{
  cf_assert(model.isNotNull());
}
      
//////////////////////////////////////////////////////////////////////
      
ICPLTE2DPuvtToCons::~ICPLTE2DPuvtToCons()
{
}

//////////////////////////////////////////////////////////////////////

void ICPLTE2DPuvtToCons::transform(const State& state, State& result)
{
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  const RealVector& refData = _icpModel->getReferencePhysicalData();
  
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal T = state[3];
  const CFreal EpR = state[4];
  const CFreal EpI = state[5];
  CFreal pdim = _icpModel->getPressureFromState(state[0])*refData[EulerTerm::P];
  CFreal Tdim = T*refData[EulerTerm::T];
  
  // set the composition
  library->setComposition(Tdim,pdim);
  library->setDensityEnthalpyEnergy(Tdim, pdim, _dhe);
  const CFreal rho = _dhe[0]/refData[EulerTerm::RHO];
  const CFreal V2 = u*u + v*v;
  
  result[0] = rho;
  result[1] = rho*u;
  result[2] = rho*v;
  result[3] = rho*(_dhe[2]/refData[EulerTerm::E] + 0.5*V2);
  result[4] = EpR;
  result[5] = EpI;
}

//////////////////////////////////////////////////////////////////////
      
void ICPLTE2DPuvtToCons::transformFromRef(const RealVector& data, State& result)
{
  // here we assume that E's are the last two components
  const CFuint firstScalarVar = _icpModel->getDataSize() - 2;
  const CFreal rho = data[EulerTerm::RHO];
  
  result[0] = rho;
  result[1] = rho*data[EulerTerm::VX];
  result[2] = rho*data[EulerTerm::VY];
  result[3] = rho*data[EulerTerm::E];
  result[4] = data[firstScalarVar];
  result[5] = data[firstScalarVar+1];
}

//////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
