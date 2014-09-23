#include "LTE.hh"
#include "Euler2DRhovtXiToPuvtLTE.hh"
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

Environment::ObjectProvider<Euler2DRhovtXiToPuvtLTE, VarSetTransformer, LTEModule, 1> 
euler2DRhovtXiToPuvtLTEProvider("Euler2DRhovtXiToPuvtLTE");

//////////////////////////////////////////////////////////////////////////////

Euler2DRhovtXiToPuvtLTE::Euler2DRhovtXiToPuvtLTE(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<EulerTerm>()),
  _dhe(3),
  _ye(),
  _xe(),
  _mmasses()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DRhovtXiToPuvtLTE::~Euler2DRhovtXiToPuvtLTE()
{
}
      
//////////////////////////////////////////////////////////////////////////////

void Euler2DRhovtXiToPuvtLTE::transform(const State& state, State& result)
{
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  const RealVector& refData = _model->getReferencePhysicalData();
  const CFuint nbSpecies = library->getNbSpecies();
  _ye.resize(nbSpecies);
  _xe.resize(nbSpecies);
  _mmasses.resize(nbSpecies);
  
  // partial densities
  for(CFuint ie = 0; ie < nbSpecies; ++ie){
    _xe[ie] = state[4 + ie];
  }
  
  library->getSpeciesMassFractions(_xe, _ye);
  library->getMolarMasses(_mmasses);
  
  CFreal sumYiOvM = 0.;
  for(CFuint ie = 0; ie < nbSpecies; ++ie) {
    sumYiOvM += _ye[ie]/_mmasses[ie];
  }
  
  const CFreal T = state[1];
  result[0] = state[0]*library->getRgas()*sumYiOvM/refData[EulerTerm::P] 
    - _model->getPressInf();
  result[1] = state[2]/refData[EulerTerm::VX];
  result[2] = state[3]/refData[EulerTerm::VY];
  result[3] = T/_model->getTempRef();
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DRhovtXiToPuvtLTE::transformFromRef(const RealVector& data, State& result)
{
  throw Common::NotImplementedException(FromHere(), "Euler2DRhovtXiToPuvtLTE::transformFromRef()");
}

//////////////////////////////////////////////////////////////////////////////

} // namespace LTE

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
