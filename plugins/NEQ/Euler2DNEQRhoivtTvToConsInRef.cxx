#include "NEQ/NEQ.hh"
#include "Euler2DNEQRhoivtTvToConsInRef.hh"
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

Environment::ObjectProvider<Euler2DNEQRhoivtTvToConsInRef, 
			    VarSetMatrixTransformer, 
			    NEQModule, 1> 
euler2DNEQRhoivtTvToConsInRefProvider("Euler2DNEQRhoivtTvToConsInRef");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQRhoivtTvToConsInRef::Euler2DNEQRhoivtTvToConsInRef
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<NEQTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQRhoivtTvToConsInRef::~Euler2DNEQRhoivtTvToConsInRef()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRhoivtTvToConsInRef::setMatrixFromRef()
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
  const CFreal u = linearData[EulerTerm::VX];
  const CFreal v = linearData[EulerTerm::VY];
  const CFuint uID  = nbSpecies;
  const CFuint vID  = nbSpecies+1;
  const CFuint eID  = nbSpecies+2;
  const CFuint evID = nbSpecies+3;
  const CFreal eT = eData->dEdT; // check this !!!!
  // DevDTv is stored in PhysicalChemicalLibrary during linearization
  const CFreal evTv = eData->dEvTv;
  const CFuint firstTv = _model->getFirstScalarVar(1);
  cf_assert(_model->getNbScalarVars(1) == 1);
  
  for (CFuint is = 0; is < nbSpecies; ++is) {
    _transMatrix(is,is) = 1.0;
    _transMatrix(uID,is) = u;
    _transMatrix(vID,is) = v;
    _transMatrix(eID,is) = linearData[EulerTerm::E];
    _transMatrix(evID,is) = linearData[firstTv];
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
