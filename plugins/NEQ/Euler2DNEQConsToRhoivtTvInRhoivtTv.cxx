#include "NEQ/NEQ.hh"
#include "Euler2DNEQConsToRhoivtTvInRhoivtTv.hh"
#include "NavierStokes/EulerPhysicalModel.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

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

Environment::ObjectProvider<Euler2DNEQConsToRhoivtTvInRhoivtTv,
			    VarSetMatrixTransformer,
			    NEQModule, 1>
euler2DNEQConsToRhoivtTvInRhoivtTvProvider("Euler2DNEQConsToRhoivtTvInRhoivtTv");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQConsToRhoivtTvInRhoivtTv::Euler2DNEQConsToRhoivtTvInRhoivtTv
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<NEQTerm>()),
  _ys(),
  _tvDim(),
  _evDim(),
  _dhe()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQConsToRhoivtTvInRhoivtTv::~Euler2DNEQConsToRhoivtTvInRhoivtTv()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQConsToRhoivtTvInRhoivtTv::setMatrix(const RealVector& state)
{
  cf_assert(_model.isNotNull());

  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();

  Common::SafePtr<PhysicalChemicalLibrary::ExtraData> eData = library->getExtraData();

  // mixture density (sum of the partial densities)
  const CFuint nbSpecies = _model->getNbScalarVars(0);
  CFreal rho = 0.0;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    rho += state[ie];
  }
  cf_assert(rho > 0.0);

  const CFreal ovRho = 1./rho;
  _ys.resize(nbSpecies);
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    _ys[ie] = state[ie]*ovRho;
  }

  // set the current species fractions in the thermodynamic library
  // this has to be done right here, before computing any other thermodynamic quantity !!!
  library->setSpeciesFractions(_ys);

  const RealVector& refData =  _model->getReferencePhysicalData();
  CFreal rhoDim = rho*refData[EulerTerm::RHO];
  CFreal T = state[nbSpecies + 2];
  CFreal Tdim = T*refData[EulerTerm::T];

  const CFuint nbTv = _model->getNbScalarVars(1);
  cf_assert(nbTv == 1);
  _dhe.resize(3 + nbTv);
  _tvDim.resize(nbTv);
  _evDim.resize(nbTv);

  const CFuint startTv = nbSpecies + 3;
  for (CFuint ie = 0; ie < nbTv; ++ie) {
    _tvDim[ie] = state[startTv + ie]*refData[EulerTerm::T];
  }

  CFreal p = library->pressure(rhoDim, Tdim, &_tvDim[0]);
  CFreal pdim = p*refData[EulerTerm::P];
  const CFreal u = state[nbSpecies];
  const CFreal v = state[nbSpecies + 1];

  // fills eData->dEdT, dEvTv, dRhoEdRhoi, dRhoEvdRhoi
  library->setDensityEnthalpyEnergy(Tdim, _tvDim, pdim,_dhe, true);

  const CFuint uID  = nbSpecies;
  const CFuint vID  = nbSpecies+1;
  const CFuint eID  = nbSpecies+2;
  const CFuint evID = nbSpecies+3;
  const CFreal eT   = eData->dEdT;
  const CFreal evTv = eData->dEvTv;
  const CFreal q = 0.5*(u*u + v*v);

  cf_assert(std::abs(eT) > 0.0);
  cf_assert(std::abs(evTv) > 0.0);

  const CFreal ovRhoET   = 1./(rho*eT);
  const CFreal ovRhoEvTv = 1./(rho*evTv);

  // Inverse of Euler2DNEQRhoivtTvToConsInRhoivtTv. Writing S = sum_j d(rho_j),
  // the forward relations invert to
  //   du   = (d(rhoU) - u S)/rho
  //   dv   = (d(rhoV) - v S)/rho
  //   dTv  = (d(rhoEv) - sum_j dRhoEvdRhoi[j] d(rho_j))/(rho evTv)
  //   dT   = (d(rhoE) - d(rhoEv) - u d(rhoU) - v d(rhoV)
  //           + sum_j (q - dRhoEdRhoi[j]) d(rho_j))/(rho eT)
  // the dRhoEvdRhoi contributions cancel in the dT row.
  _transMatrix = 0.0;

  for (CFuint js = 0; js < nbSpecies; ++js) {
    _transMatrix(js,js)   = 1.0;
    _transMatrix(uID,js)  = -u*ovRho;
    _transMatrix(vID,js)  = -v*ovRho;
    _transMatrix(eID,js)  = (q - eData->dRhoEdRhoi[js])*ovRhoET;
    _transMatrix(evID,js) = -eData->dRhoEvdRhoi[js]*ovRhoEvTv;
  }

  _transMatrix(uID,uID) = ovRho;
  _transMatrix(vID,vID) = ovRho;

  _transMatrix(eID,uID)  = -u*ovRhoET;
  _transMatrix(eID,vID)  = -v*ovRhoET;
  _transMatrix(eID,eID)  =  ovRhoET;
  _transMatrix(eID,evID) = -ovRhoET;

  _transMatrix(evID,evID) = ovRhoEvTv;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
