#include "NEQ/NEQ.hh"
#include "Euler2DNEQRhoivtTvToConsInRhoivtTv.hh"
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

Environment::ObjectProvider<Euler2DNEQRhoivtTvToConsInRhoivtTv, 
			    VarSetMatrixTransformer, 
			    NEQModule, 1> 
euler2DNEQRhoivtTvToConsInRhoivtTvProvider("Euler2DNEQRhoivtTvToConsInRhoivtTv");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQRhoivtTvToConsInRhoivtTv::Euler2DNEQRhoivtTvToConsInRhoivtTv
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

Euler2DNEQRhoivtTvToConsInRhoivtTv::~Euler2DNEQRhoivtTvToConsInRhoivtTv()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRhoivtTvToConsInRhoivtTv::setMatrix(const RealVector& state)
{
  cf_assert(_model.isNotNull());
 
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  Common::SafePtr<PhysicalChemicalLibrary::ExtraData> eData = library->getExtraData();
     
  // Set the mixture density (sum of the partial densities)
  const CFuint nbSpecies = _model->getNbScalarVars(0);
  CFreal rho = 0.0;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    rho += state[ie];
  }
  
  // Set the species
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
  
  // set the vibrational temperature
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
  
  // unused //  const CFreal ovHref = 1./refData[EulerTerm::H];
  // vibrational temperatures
  library->setDensityEnthalpyEnergy(Tdim, _tvDim, pdim,_dhe, true);
  
  const CFuint uID  = nbSpecies;
  const CFuint vID  = nbSpecies+1;
  const CFuint eID  = nbSpecies+2;
  const CFuint evID = nbSpecies+3;
  const CFreal eT = eData->dEdT; // check this
  // DevDTv is stored in PhysicalChemicalLibrary during linearization
  const CFreal evTv = eData->dEvTv;
  const CFreal q = 0.5*(u*u + v*v);		
  
  for (CFuint is = 0; is < nbSpecies; ++is) {
    _transMatrix(is,is) = 1.0;
    _transMatrix(uID,is) = u;
    _transMatrix(vID,is) = v;
    _transMatrix(eID,is) = eData->dRhoEdRhoi[is] + eData->dRhoEvdRhoi[is] + q;
    _transMatrix(evID,is) = eData->dRhoEvdRhoi[is];
  }
  
  _transMatrix(uID,uID) = rho;
  _transMatrix(vID,vID) = rho;
  
  _transMatrix(eID,uID)  = rho*u;
  _transMatrix(eID,vID)  = rho*v;
  _transMatrix(eID,eID)  = rho*eT;
  _transMatrix(eID,evID) = rho*evTv;
  
  _transMatrix(evID,evID) = rho*evTv;
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
