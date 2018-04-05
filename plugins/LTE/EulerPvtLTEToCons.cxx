#include "LTE.hh"
#include "EulerPvtLTEToCons.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "NavierStokes/EulerTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace LTE {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<EulerPvtLTEToCons, VarSetTransformer, LTEModule, 1> 
euler1DPvtLTEToConsProvider("Euler1DPvtLTEToCons");

Environment::ObjectProvider<EulerPvtLTEToCons, VarSetTransformer, LTEModule, 1> 
euler2DPuvtLTEToConsProvider("Euler2DPuvtLTEToCons");

Environment::ObjectProvider<EulerPvtLTEToCons, VarSetTransformer, LTEModule, 1> 
euler3DPvtLTEToConsProvider("Euler3DPvtLTEToCons");

// 2D and 1/2
Environment::ObjectProvider<EulerPvtLTEToCons, VarSetTransformer, LTEModule, 1> 
euler2DPvtLTEToConsProvider("Euler2DPvtLTEToCons");

//////////////////////////////////////////////////////////////////////////////

EulerPvtLTEToCons::EulerPvtLTEToCons(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>()), 
  m_rho(0.),
  _dhe(3)
{
}

//////////////////////////////////////////////////////////////////////////////

EulerPvtLTEToCons::~EulerPvtLTEToCons()
{
}
      
//////////////////////////////////////////////////////////////////////////////

void EulerPvtLTEToCons::transform(const State& state, State& result)
{
  const RealVector& refData = _model->getReferencePhysicalData();
  
  Common::SafePtr<PhysicalModelImpl> pm =
    PhysicalModelStack::getActive()->getImplementor();
  const CFuint dim = (pm->is2DHalf()) ?
    3 : PhysicalModelStack::getActive()->getDim();
  const CFuint TID = dim+1;
  
  CFreal p = _model->getPressureFromState(state[0]);
  CFreal T = state[TID];
  CFreal pdim = p*refData[EulerTerm::P];
  CFreal Tdim = T*_model->getTempRef();
  
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    pm->getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  // set the composition
  library->setComposition(Tdim,pdim);
  library->setDensityEnthalpyEnergy(Tdim, pdim,_dhe);
  
  m_rho = _dhe[0]/refData[EulerTerm::RHO];
  result[0] = m_rho;
  
  CFreal V2 = 0.;
  for (CFuint i=0; i < dim; ++i) {
    const CFuint uID = i+1;
    const CFreal u = state[uID];
    result[uID] = m_rho*u; 
    V2 += u*u;
  }
  
  result[TID] = m_rho*(_dhe[2]/refData[EulerTerm::H] + 0.5*V2);
}
      
//////////////////////////////////////////////////////////////////////////////

void EulerPvtLTEToCons::transformFromRef(const RealVector& data, State& result)
{
  m_rho = data[EulerTerm::RHO];
  result[0] = m_rho;
  
  Common::SafePtr<PhysicalModelImpl> pm =
    PhysicalModelStack::getActive()->getImplementor();
  const CFuint dim = (pm->is2DHalf()) ?
    3 : PhysicalModelStack::getActive()->getDim();
  for (CFuint i=0; i < dim; ++i) {
    result[i+1] = m_rho*data[EulerTerm::VX+i];
  }
  
  const CFuint TID = dim+1;
  result[TID] = m_rho*data[EulerTerm::E];
}

//////////////////////////////////////////////////////////////////////////////

} // namespace LTE

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
