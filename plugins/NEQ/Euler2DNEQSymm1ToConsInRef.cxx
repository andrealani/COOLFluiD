#include "NEQ/NEQ.hh"
#include "Euler2DNEQSymm1ToConsInRef.hh"
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

Environment::ObjectProvider<Euler2DNEQSymm1ToConsInRef, 
			    VarSetMatrixTransformer, 
			    NEQModule, 1> 
euler2DNEQSymm1ToConsInRefProvider("Euler2DNEQSymm1ToConsInRef");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQSymm1ToConsInRef::Euler2DNEQSymm1ToConsInRef
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<NEQTerm>()), _ys(), _alpha(), _RiGas(), _mmasses(), _fcoeff()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQSymm1ToConsInRef::~Euler2DNEQSymm1ToConsInRef()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm1ToConsInRef::setMatrixFromRef()
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
  // const CFreal ovrho = 1./rho;
  const CFreal u = linearData[EulerTerm::VX];
  const CFreal v = linearData[EulerTerm::VY];
  const CFreal H = linearData[EulerTerm::H];
  const CFreal a = linearData[EulerTerm::A];
  const CFreal a2 = linearData[EulerTerm::A]*linearData[EulerTerm::A];
  const CFreal ova = 1./a;
  const CFreal ova2 = 1./a2;
  const CFuint uID  = nbSpecies;
  const CFuint vID  = nbSpecies+1;
  const CFuint eID  = nbSpecies+2;
  const CFuint evID = nbSpecies+3;
  const CFuint firstTv = _model->getFirstScalarVar(1);
  cf_assert(_model->getNbScalarVars(1) == 1);
  
  const CFuint firstSpecies = _model->getFirstScalarVar(0);
  CFreal sumAlphaYs = 0.0;

  _ys.resize(nbSpecies);
  _alpha.resize(nbSpecies);
  _RiGas.resize(nbSpecies);
  _mmasses.resize(nbSpecies);
  _fcoeff.resize(nbSpecies);
  
  //The previous vectors are initialized:

  cf_assert(_ys.size() == nbSpecies);
  cf_assert(_alpha.size() == nbSpecies);
  //  cf_assert(p > 0.0);
  //cf_assert(T > 0.0);
  //cf_assert(rho > 0.0);
  
  library->setRiGas(_RiGas);
  library->getMolarMasses(_mmasses);
  
  vector<CFuint> moleculeIDs;
  library->setMoleculesIDs(moleculeIDs);
  vector<bool> flag(nbSpecies, false);
  for (CFuint i = 0; i < moleculeIDs.size(); ++i) {
    flag[moleculeIDs[i]] = true;
  }
  
  for (CFuint i = 0; i < nbSpecies; ++i) {
    _fcoeff[i] = (flag[i]) ? 2.5 : 1.5;
  }
  
  CFreal numBeta = 0.;
  CFreal denBeta = 0.;
  const CFuint start = (library->presenceElectron()) ? 1 : 0;
  for (CFuint i = start; i < nbSpecies; ++i) {
    const CFreal sigmai = linearData[firstSpecies + i]/_mmasses[i];
    numBeta += sigmai;
    denBeta += sigmai*_fcoeff[i];
  }
  
  const CFreal T = linearData[EulerTerm::T];
  const CFreal beta = numBeta/denBeta;   
   
  for (CFuint is = 0; is < nbSpecies; ++is) {
    _ys[is] = linearData[firstSpecies + is];
    if (_ys[is] > 1.1) {
      cout << "_ys > 1.1 = " << _ys << endl;
      // cf_assert(_ys[is] <= 1.1);
    }
  }
  
  CFreal phi = 0.0; 
  if (library->presenceElectron()) {
    // assume that electrons have 0 as mixture ID  
    phi = _RiGas[0]*_ys[0]/eData->dEvTv - beta;
  }
  else {
    phi =-beta;
  }
  
  const CFreal V2=linearData[EulerTerm::V]*linearData[EulerTerm::V];
  const CFreal q = 0.5*V2;
  const CFreal bq = beta*q;
  CFreal Tq = 0.0;
  
  for (CFuint is = 0; is < nbSpecies; ++is) {
    // here T must be substituted with Tv if there is ionization
    // here phi has to be changed for the ionized case
    if (!library->presenceElectron()) {
      Tq = T;
    }
    else {
      if (is == 0) {
	Tq = 0.0; //To implement: get vibrational temperature
      }
      else {
	Tq = T;
      }
    }
    
    _alpha[is] = _RiGas[is]*Tq + bq - beta*((eData->energyTr)[is] + (eData->energyVib)[is]) - phi*((eData->energyVib)[is]);
    /*cout<<"Index is = " << is <<endl;
    cout<<"Tq = " << Tq <<endl;
    cout<<"Tq - T = " << Tq-T <<endl;
    cout<<"_alpha[is] = " << _alpha[is]<<endl;
    cout<<"_RiGas[is]*Tq + bq - beta*(eData->energyTr)[is] = " << _RiGas[is]*Tq + bq - beta*(eData->energyTr)[is] <<endl;
    cout.precision(20);
    cout<<"difference = " << _RiGas[is]*Tq + bq - beta*(eData->energyTr)[is] - _alpha[is] <<endl;*/
    //The previous commented portion of code shows a difference of the order of  1e-10.
    //Thus we comment out the next assertion:
    //cf_assert(_RiGas[is]*Tq + bq - beta*(eData->energyTr)[is] == _alpha[is]);
    sumAlphaYs += _alpha[is]*_ys[is];
  }
  
  // Vectors initialization ends here.
  const CFreal ev = linearData[firstTv];
  
  for (CFuint is = 0; is < nbSpecies; ++is) {
    _transMatrix(is,is) = -ova2;
    _transMatrix(uID,is) = -ova2*u;
    _transMatrix(vID,is) = -ova2*v;
    _transMatrix(eID,is) = ova2*(phi*ev/beta-V2)+ova2*_alpha[is]/beta;
    _transMatrix(evID,is) = -ova2*ev;
  }
  
  for (CFuint is = 0; is < nbSpecies; ++is) {
    _transMatrix(is,eID) = ova*rho*_ys[is];
  }
  
  _transMatrix(uID,uID) = rho;
  _transMatrix(uID,vID) = 0.0;
  _transMatrix(uID,eID) = rho*ova*u;
  _transMatrix(uID,evID) = 0.0;
  
  _transMatrix(vID,uID) = 0.0;
  _transMatrix(vID,vID) = rho;
  _transMatrix(vID,eID) = rho*ova*v;
  _transMatrix(vID,evID) = 0.0;

  _transMatrix(eID,uID)  = rho*u;
  _transMatrix(eID,vID)  = rho*v;
  _transMatrix(eID,eID)  = rho*ova*(H - phi*ev/beta);//ova*rho*(V2-2*phi*ev/beta)+rho*a/beta-ova*rho*sumAlphaYs/beta;
  _transMatrix(eID,evID) = ova2*rho*phi/beta;
  
  _transMatrix(evID,uID) = 0.0;
  _transMatrix(evID,vID) = 0.0;
  _transMatrix(evID,eID) = 2*rho*ev*ova;
  _transMatrix(evID,evID) = -rho*ova2;
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
