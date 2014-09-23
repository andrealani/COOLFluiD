#include "NEQ/NEQ.hh"
#include "Euler2DNEQConsToSymmInRef.hh"
#include "NavierStokes/EulerPhysicalModel.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
//
#include "MathTools/LUInverter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DNEQConsToSymmInRef, 
			    VarSetMatrixTransformer, 
			    NEQModule, 1> 
euler2DNEQConsToSymmInRefProvider("Euler2DNEQConsToSymmInRef");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQConsToSymmInRef::Euler2DNEQConsToSymmInRef
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<NEQTerm>()), 
  _ys(), _alpha(), _RiGas(), _mmasses()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQConsToSymmInRef::~Euler2DNEQConsToSymmInRef()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQConsToSymmInRef::setMatrixFromRef()
{
//  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
//   RealMatrix direct(nbEqs, nbEqs);
//   RealMatrix inverse(nbEqs, nbEqs);
  
//   computeDirect(direct);
//   computeInverse(inverse);
  
//   cout << "INVERSE 1"<< endl << inverse << endl;
//   RealMatrix inverseMat(nbEqs, nbEqs);
//   LUInverter lu(nbEqs);
//   lu.invert(direct, inverseMat);
  
//   cout << "INVERSE 2"<< endl <<inverseMat <<  endl;
  
  
//   RealMatrix product(nbEqs, nbEqs);
//   product = direct*inverseMat;
//   cout << "D*I" << endl << product << endl;
  
//   product = inverseMat*direct;
//   cout << "I*D" << endl << product << endl;
//   abort();
  
  cf_assert(_model.isNotNull());
  const RealVector& linearData = _model->getPhysicalData();
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  Common::SafePtr<PhysicalChemicalLibrary::ExtraData> eData = library->getExtraData();
  
  const CFuint nbSpecies = _model->getNbScalarVars(0);
  const CFreal rho = linearData[EulerTerm::RHO];
  const CFreal ovrho = 1./rho;
  const CFreal u = linearData[EulerTerm::VX];
  const CFreal v = linearData[EulerTerm::VY];
  const CFreal a = linearData[EulerTerm::A];
  const CFreal a2 = linearData[EulerTerm::A]*linearData[EulerTerm::A];
  const CFreal ova = 1./a;
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
  
  cf_assert(_ys.size() == nbSpecies);
  cf_assert(_alpha.size() == nbSpecies);
  
  library->setRiGas(_RiGas);
  library->getMolarMasses(_mmasses);
  SafePtr<RealVector> fcoeff = library->getAtomicityCoeff();
  
  CFreal numBeta = 0.;
  CFreal denBeta = 0.;
  const CFuint start = (library->presenceElectron()) ? 1 : 0;
  for (CFuint i = start; i < nbSpecies; ++i) {
    const CFreal sigmai = linearData[firstSpecies + i]/_mmasses[i];
    numBeta += sigmai;
    denBeta += sigmai*(*fcoeff)[i];
  }
  
  const CFreal T = linearData[EulerTerm::T];
  const CFreal beta = numBeta/denBeta;
  //const CFreal p = linearData[EulerTerm::P];
  //const CFreal cvTr = eData->dEdT;
  //const CFreal beta = p/(rho*T)/cvTr;
   
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
    for (CFuint js = 0; js < nbSpecies; ++js) {
      	_transMatrix(is,js) = _ys[is]*_alpha[js];				
      if (is==js){
	_transMatrix(is,js)-=a2;
      }     
    }

    _transMatrix(uID,is) = -ovrho*u;
    _transMatrix(vID,is) = -ovrho*v;
    _transMatrix(eID,is) = ovrho*ova*_alpha[is];
    _transMatrix(evID,is) = ovrho*ev*_alpha[is];
  }
  
    for (CFuint is = 0; is < nbSpecies; ++is) {
    _transMatrix(is,uID) = -beta*_ys[is]*u;
    _transMatrix(is,vID) = -beta*_ys[is]*v;
    _transMatrix(is,eID) = beta*_ys[is];
    _transMatrix(is,evID) = phi*_ys[is];
  }

  _transMatrix(uID,uID) = ovrho;
  _transMatrix(uID,vID) = 0.0;
  _transMatrix(uID,eID) = 0.0;
  _transMatrix(uID,evID) = 0.0;

  _transMatrix(vID,uID) = 0.0;
  _transMatrix(vID,vID) = ovrho;
  _transMatrix(vID,eID) = 0.0;
  _transMatrix(vID,evID) = 0.0;
  
  _transMatrix(eID,uID) = -ovrho*ova*beta*u;
  _transMatrix(eID,vID) = -ovrho*ova*beta*v;
  _transMatrix(eID,eID)  = ovrho*ova*beta;
  _transMatrix(eID,evID) = ovrho*ova*phi;
  
  _transMatrix(evID,uID) = -ovrho*beta*ev*u;
  _transMatrix(evID,vID) = -ovrho*beta*ev*v;
  _transMatrix(evID,eID) = ovrho*beta*ev;      
  _transMatrix(evID,evID) = ovrho*(phi*ev-a2);//_transMatrix checked on 2009/3/30.
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQConsToSymmInRef::computeDirect(RealMatrix& mat)
{  
  //Note: this method contains the definition of the matrix consToSymm:
  cf_assert(_model.isNotNull());
  const RealVector& linearData = _model->getPhysicalData();
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  Common::SafePtr<PhysicalChemicalLibrary::ExtraData> eData = library->getExtraData();
  
  const CFuint nbSpecies = _model->getNbScalarVars(0);
  const CFreal rho = linearData[EulerTerm::RHO];
  const CFreal ovrho = 1./rho;
  const CFreal u = linearData[EulerTerm::VX];
  const CFreal v = linearData[EulerTerm::VY];
  const CFreal a = linearData[EulerTerm::A];
  const CFreal a2 = linearData[EulerTerm::A]*linearData[EulerTerm::A];
  const CFreal ova = 1./a;
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

  cf_assert(_ys.size() == nbSpecies);
  cf_assert(_alpha.size() == nbSpecies);
  
  library->setRiGas(_RiGas);
  library->getMolarMasses(_mmasses);
  SafePtr<RealVector> fcoeff = library->getAtomicityCoeff();
  
  CFreal numBeta = 0.;
  CFreal denBeta = 0.;
  const CFuint start = (library->presenceElectron()) ? 1 : 0;
  for (CFuint i = start; i < nbSpecies; ++i) {
    const CFreal sigmai = linearData[firstSpecies + i]/_mmasses[i];
    numBeta += sigmai;
    denBeta += sigmai*(*fcoeff)[i];
  }
  
  const CFreal T = linearData[EulerTerm::T];
  const CFreal beta = numBeta/denBeta;
  //const CFreal p = linearData[EulerTerm::P];
  //const CFreal cvTr = eData->dEdT;
  //const CFreal beta = p/(rho*T)/cvTr;
   
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
    for (CFuint js = 0; js < nbSpecies; ++js) {
      	mat(is,js) = _ys[is]*_alpha[js];				
      if (is==js){
	mat(is,js)-=a2;
      }     
    }

    mat(uID,is) = -ovrho*u;
    mat(vID,is) = -ovrho*v;
    mat(eID,is) = ovrho*ova*_alpha[is];
    mat(evID,is) = ovrho*ev*_alpha[is];
  }
  
    for (CFuint is = 0; is < nbSpecies; ++is) {
    mat(is,uID) = -beta*_ys[is]*u;
    mat(is,vID) = -beta*_ys[is]*v;
    mat(is,eID) = beta*_ys[is];
    mat(is,evID) = phi*_ys[is];
  }

  mat(uID,uID) = ovrho;
  mat(uID,vID) = 0.0;
  mat(uID,eID) = 0.0;
  mat(uID,evID) = 0.0;

  mat(vID,uID) = 0.0;
  mat(vID,vID) = ovrho;
  mat(vID,eID) = 0.0;
  mat(vID,evID) = 0.0;
  
  mat(eID,uID) = -ovrho*ova*beta*u;
  mat(eID,vID) = -ovrho*ova*beta*v;
  mat(eID,eID)  = ovrho*ova*beta;
  mat(eID,evID) = ovrho*ova*phi;
  
  mat(evID,uID) = -ovrho*beta*ev*u;
  mat(evID,vID) = -ovrho*beta*ev*v;
  mat(evID,eID) = ovrho*beta*ev;      
  mat(evID,evID) = ovrho*(phi*ev-a2);//mat checked on 2009/3/30.
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQConsToSymmInRef::computeInverse(RealMatrix& mat)
{
  //Note: this method computes the analytical symmToCons
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
  
  cf_assert(_ys.size() == nbSpecies);
  cf_assert(_alpha.size() == nbSpecies);
  
  library->setRiGas(_RiGas);
  library->getMolarMasses(_mmasses);
  SafePtr<RealVector> fcoeff = library->getAtomicityCoeff();
 
  CFreal numBeta = 0.;
  CFreal denBeta = 0.;
  const CFuint start = (library->presenceElectron()) ? 1 : 0;
  for (CFuint i = start; i < nbSpecies; ++i) {
    const CFreal sigmai = linearData[firstSpecies + i]/_mmasses[i];
    numBeta += sigmai;
    denBeta += sigmai*(*fcoeff)[i];
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
    mat(is,is) = -ova2;
    mat(uID,is) = -ova2*u;
    mat(vID,is) = -ova2*v;
    mat(eID,is) = ova2*_alpha[is]/beta-ova2*V2;
    mat(evID,is) = 0.0;
  }
  
for (CFuint is = 0; is < nbSpecies; ++is) {
  mat(is,eID) = ova*rho*_ys[is];
 }

  mat(uID,uID) = rho;
  mat(uID,vID) = 0.0;
  mat(uID,eID) = ova*rho*u;
  mat(uID,evID) = 0.0;
  
  mat(vID,uID) = 0.0;
  mat(vID,vID) = rho;
  mat(vID,eID) = ova*rho*v;
  mat(vID,evID) = 0.0;

  mat(eID,uID)  = rho*u;
  mat(eID,vID)  = rho*v;
  mat(eID,eID)  = ova*rho*H;
  mat(eID,evID) = ova2*phi*rho/beta;
  
  mat(evID,uID) = 0.0;
  mat(evID,vID) = 0.0;
  mat(evID,eID) = ova*rho*ev;
  mat(evID,evID) = -ova2*rho;//mat checked on 2009/4/2. Error in position (eID,eID) fixed.
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
