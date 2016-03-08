#include "NEQ/NEQ.hh"
#include "NEQ/Euler2DNEQRhoivtTvToPivtTv.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DNEQRhoivtTvToPivtTv, VarSetTransformer, NEQModule, 1>
euler2DNEQRhoivtTvToPivtTvProvider("Euler2DNEQRhoivtTvToPivtTv");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQRhoivtTvToPivtTv::Euler2DNEQRhoivtTvToPivtTv
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<NEQTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQRhoivtTvToPivtTv::~Euler2DNEQRhoivtTvToPivtTv()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRhoivtTvToPivtTv::transform(const State& state, State& result)
{
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  const CFuint nbSpecies = _model->getNbScalarVars(0);
  RealVector ri(nbSpecies);
  library->setRiGas(ri);
  
  const CFreal T = state[nbSpecies + 2]; 
  // partial pressures for heavy particles
  for(CFuint ie = 1; ie < nbSpecies; ++ie) {
    result[ie] = state[ie]*ri[ie]*T;
    if (result[ie] < -1e-14) {
      CFLog(DEBUG_MIN, "density[" << ie  << "] = " << state[ie] << " => pressure[" << ie  << "] = " << result[ie] << "\n");
    }
  }
  
  // electron temperature is treated separately 
  // AL: we assume electrons to be the first!
  const CFuint nbTv = _model->getNbScalarVars(1);
  const CFuint nbTe = library->getNbTe();
  const CFuint nbTvH = nbTv - nbTe; 
  const CFuint startTv = nbSpecies + 3;
  
  // partial pressure for electrons
  if (nbTe == 1) {
    cf_assert(startTv + nbTvH < state.size());
    result[0] = state[0]*ri[0]*state[startTv + nbTvH];
  }
  else {
    // AL: we assume that the first Tv is associated to electrons
    cf_assert(startTv < state.size());
    result[0] = state[0]*ri[0]*state[startTv];
  } 
  if (result[0] < -1e-14) {
    CFLog(DEBUG_MIN, "density[0] = " << state[0] << " => pressure[0] = " << result[0] << "\n");
  }
  
  const CFuint nbEqs = nbSpecies+3+nbTv; 
  for(CFuint ie = nbSpecies; ie < nbEqs; ++ie) {
    result[ie] = state[ie];
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRhoivtTvToPivtTv::transformFromRef(const RealVector& data, State& result)
{
  throw Common::NotImplementedException(FromHere(), "Euler2DNEQRhoivtTvToPivtTv::transformFromRef()");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
