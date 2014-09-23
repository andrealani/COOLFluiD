#include "NEQ/NEQ.hh"
#include "Euler2DNEQRoeToConsInRef.hh"
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

Environment::ObjectProvider<Euler2DNEQRoeToConsInRef, VarSetMatrixTransformer, 
NEQModule, 1> euler2DNEQRoeToConsInRefProvider("Euler2DNEQRoeToConsInRef");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQRoeToConsInRef::Euler2DNEQRoeToConsInRef
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<NEQTerm>()),
  _RiGas()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQRoeToConsInRef::~Euler2DNEQRoeToConsInRef()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRoeToConsInRef::setMatrixFromRef()
{
  cf_assert(_model.isNotNull());
  const RealVector& linearData = _model->getPhysicalData();
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  Common::SafePtr<PhysicalChemicalLibrary::ExtraData> eData = library->getExtraData();
  
  const CFuint nbSpecies = _model->getNbScalarVars(0);
  const CFreal rho = linearData[EulerTerm::RHO];
  const CFreal sqRho = sqrt(rho);
  const CFreal p = linearData[EulerTerm::P];
  const CFreal T = linearData[EulerTerm::T];
  const CFreal H = linearData[EulerTerm::H];
  const CFreal u = linearData[EulerTerm::VX];
  const CFreal v = linearData[EulerTerm::VY];
  const CFreal V2 = linearData[EulerTerm::V]*linearData[EulerTerm::V];
  const CFuint uID  = nbSpecies;
  const CFuint vID  = nbSpecies+1;
  const CFuint eID  = nbSpecies+2;
  const CFuint evID = nbSpecies+3;
  const CFreal cvTr = eData->dEdT; // check this !!!!
  const CFuint firstTv = _model->getFirstScalarVar(1);
  const CFuint firstSpecies = _model->getFirstScalarVar(0);
  const CFreal ev = linearData[firstTv];
  const CFreal dd = rho*T*cvTr + p;
  const CFreal Rb = p/(sqRho*T);
  const CFreal c = sqRho*cvTr + Rb;
  const CFreal c2 = c*c;
  const CFreal dRhoEdzev = sqRho*p/dd;

  _RiGas.resize(nbSpecies);
  library->setRiGas(_RiGas); 
  
  for (CFuint is = 0; is < nbSpecies; ++is) {
    _transMatrix(is,is) = sqRho*(linearData[firstSpecies + is] + 1.);
    _transMatrix(uID,is) = sqRho*u;
    _transMatrix(vID,is) = sqRho*v;
    
    const CFreal dcDZi = cvTr + _RiGas[is];
    const CFreal cH  = H*cvTr*sqRho*(2*c - sqRho*dcDZi);
    const CFreal cEv = ev*((Rb - _RiGas[is]*sqRho)*c - Rb*sqRho*dcDZi);
    const CFreal cV  = 0.5*V2*sqRho*(_RiGas[is]*c - Rb*dcDZi);
    _transMatrix(eID,is) = sqRho/c2*(cH + cEv + cV);
    _transMatrix(evID,is) = sqRho*ev;
  }
  
  _transMatrix(uID,uID) = sqRho;
  _transMatrix(vID,vID) = sqRho;
  
  _transMatrix(eID,uID)  = dRhoEdzev*u;
  _transMatrix(eID,vID)  = dRhoEdzev*v;
  _transMatrix(eID,eID)  = (rho*sqRho*cvTr*T)/dd;
  _transMatrix(eID,evID) = dRhoEdzev;
  
  _transMatrix(evID,evID) = sqRho;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
