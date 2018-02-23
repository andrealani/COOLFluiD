#include "ICP/ICP.hh"
#include "ICP/ICPLTEPvtToCons.hh"
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

Environment::ObjectProvider<ICPLTEPvtToCons, VarSetTransformer, ICPModule,1>
icpLTE2DPuvtToConsProvider("ICPLTE2DPuvtToCons");

Environment::ObjectProvider<ICPLTEPvtToCons, VarSetTransformer, ICPModule,1>
icpLTE2DPvtToConsProvider("ICPLTE2DPvtToCons");
      
//////////////////////////////////////////////////////////////////////

ICPLTEPvtToCons::ICPLTEPvtToCons
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _icpModel(model->getConvectiveTerm().d_castTo<PTERM>()),
  _dhe(3)
{
  cf_assert(model.isNotNull());
}
      
//////////////////////////////////////////////////////////////////////
      
ICPLTEPvtToCons::~ICPLTEPvtToCons()
{
}

//////////////////////////////////////////////////////////////////////

void ICPLTEPvtToCons::transform(const State& state, State& result)
{
  Common::SafePtr<PhysicalModelImpl> pm =
    PhysicalModelStack::getActive()->getImplementor();
  const CFuint dim = (pm->is2DHalf()) ?
    3 : PhysicalModelStack::getActive()->getDim();
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    pm->getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  const RealVector& refData = _icpModel->getReferencePhysicalData();
  
  const CFuint TID = dim+1;
  const CFuint erID = TID+1;
  const CFuint eiID = TID+2;
  const CFreal T = state[TID];
  const CFreal EpR = state[erID];
  const CFreal EpI = state[eiID];
  CFreal pdim = _icpModel->getPressureFromState(state[0])*refData[EulerTerm::P];
  CFreal Tdim = T*refData[EulerTerm::T];
  
  // set the composition
  library->setComposition(Tdim,pdim);
  library->setDensityEnthalpyEnergy(Tdim, pdim, _dhe);
  const CFreal rho = _dhe[0]/refData[EulerTerm::RHO];
  result[0] = rho;

  CFreal V2 = 0.;
  for (CFuint i = 1; i <= dim; ++i) {
    const CFreal ui = state[i];
    result[i] = rho*ui;
    V2 += ui*ui;
  }
  
  result[TID] = rho*(_dhe[2]/refData[EulerTerm::E] + 0.5*V2);
  result[erID] = EpR;
  result[eiID] = EpI;
}

//////////////////////////////////////////////////////////////////////
      
void ICPLTEPvtToCons::transformFromRef(const RealVector& data, State& result)
{
  Common::SafePtr<PhysicalModelImpl> pm =
    PhysicalModelStack::getActive()->getImplementor();
  const CFuint dim = (pm->is2DHalf()) ?
    3 : PhysicalModelStack::getActive()->getDim();
  
  // here we assume that E's are the last two components
  const CFuint firstScalarVar = _icpModel->getDataSize() - 2;
  const CFreal rho = data[EulerTerm::RHO];
  result[0] = rho;
  for (CFuint i = 0; i < dim; ++i) {
    result[i+1] = rho*data[EulerTerm::VX+i];
  }
  result[dim+1] = rho*data[EulerTerm::E];
  result[dim+2] = data[firstScalarVar];
  result[dim+3] = data[firstScalarVar+1];
}

//////////////////////////////////////////////////////////////////////

    } // namespace ICP

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
