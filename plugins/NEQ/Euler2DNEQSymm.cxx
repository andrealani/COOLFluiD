#include "NEQ.hh"
#include "Euler2DNEQSymm.hh"
#include "Common/NotImplementedException.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DNEQSymm, ConvectiveVarSet, NEQModule, 1>
euler2DNEQSymmProvider("Euler2DNEQSymm");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQSymm::Euler2DNEQSymm(Common::SafePtr<BaseTerm> term) :
  MultiScalarVarSet<Euler2DVarSet>(term),//This sequence is an initialization.
  _library(CFNULL),
  _Rgas(),
  _dhe(3),
  _ys(),
  _rightEv(),
  _leftEv(),
  _alpha(),
  _RiGas(),
  _mmasses(),
  _fcoeff()
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint nbTv = getModel()->getNbScalarVars(1);
  
  vector<std::string> names(nbSpecies + 3 + nbTv);
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    names[ie] = "rho" + StringOps::to_str(ie);
  }
  names[nbSpecies]     = "rhoU";
  names[nbSpecies + 1] = "rhoV";
  names[nbSpecies + 2] = "rhoE";
  
  const CFuint startTv = nbSpecies + 3; 
  for (CFuint ie = 0; ie < nbTv; ++ie) {
    names[startTv + ie] = "rhoEv" + StringOps::to_str(ie);
  }
  
  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQSymm::~Euler2DNEQSymm()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm::computeEigenValuesVectors(RealMatrix& rightEv,
                                           RealMatrix& leftEv,
                                           RealVector& eValues,
                                           const RealVector& normal)
{
  throw Common::NotImplementedException(FromHere(), "Euler2DNEQSymm::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler2DNEQSymm::getBlockSeparator() const
{
  return getModel()->getNbScalarVars(0);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm::splitJacobian(RealMatrix& jacobPlus,
				   RealMatrix& jacobMin,
				   RealVector& eValues,
				   const RealVector& normal)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint nbEqs = nbSpecies+4;//nbSp + u + v + p + ev
  const CFuint jacobSize = jacobPlus.nbRows();
  
  if (jacobSize == nbEqs-nbSpecies-1){
    splitJacobian_FullDecoupling(jacobPlus,
				 jacobMin,
				 eValues,
				 normal);
  }
  
  else if (jacobSize == nbEqs-nbSpecies){
    splitJacobian_PartialDecoupling(jacobPlus,
				    jacobMin,
				    eValues,
				    normal);
  }

  else if (jacobSize == nbEqs) {
    splitJacobian_NoDecoupling(jacobPlus,
			       jacobMin,
			       eValues,
			       normal);
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm::computePhysicalData(const State& state,
					 RealVector& data)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DNEQSymm::computePhysicalData()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm::computeStateFromPhysicalData(const RealVector& data,
					 State& state)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DNEQSymm::computeStateFromPhysicalData()");
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler2DNEQSymm::getSpeed(const State& state) const
{
  const CFreal u = state[1]/state[0];
  const CFreal v = state[2]/state[0];
  return sqrt(u*u + v*v);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm::setDimensionalValues(const State& state,
                                          RealVector& result)
{
  throw Common::NotImplementedException
      (FromHere(), "Euler2DNEQSymm::setDimensionalValues()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm::setAdimensionalValues(const Framework::State& state,
                             RealVector& result)
{
  throw Common::NotImplementedException
      (FromHere(), "Euler2DNEQSymm::setAdimensionalValues() not implemented");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm::setDimensionalValuesPlusExtraValues
(const State& state, RealVector& result,
 RealVector& extra)
{
  throw Common::NotImplementedException
    (FromHere(), "Euler2DNEQSymm::setDimensionalValuesPlusExtraValues()");
}
      
//////////////////////////////////////////////////////////////////////////////

vector<std::string> Euler2DNEQSymm::getExtraVarNames() const
{
  cf_assert(_library.isNotNull());

  vector<std::string> names(3);
  names[0] = "rho";
  names[1] = "H";
  names[2] = "M";

  return names;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm::setup()
{
  MultiScalarVarSet<Euler2DVarSet>::setup();
  
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  // set the equation set data for each of the equation subsets
  // first equation subset
  MultiScalarVarSet<Euler2DVarSet>::getEqSetData().resize(2);
  MultiScalarVarSet<Euler2DVarSet>::getEqSetData()[0].setup(0,0,nbSpecies);
  
  // second equation subset
  Euler2DVarSet::getEqSetData().resize(1);
  Euler2DVarSet::getEqSetData()[0].setup(1,nbSpecies,3);
  
  const CFuint nbTv = getModel()->getNbScalarVars(1);
  // third equation subset
  MultiScalarVarSet<Euler2DVarSet>::getEqSetData()[1].setup(2,nbSpecies + 3,nbTv);
  
  _library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(_library.isNotNull());
  
  _Rgas = _library->getRgas();
  
  _dhe.resize(3 + nbTv);
  _ys.resize(nbSpecies);
  const CFuint totNbEqs = 3 + nbTv + nbSpecies;
   
  _rightEv.resize(totNbEqs, totNbEqs);
  _rightEv = 0.0;
  
  _leftEv.resize(totNbEqs, totNbEqs);
  _leftEv = 0.0;
  
  _alpha.resize(nbSpecies);
  _RiGas.resize(nbSpecies);
  _library->setRiGas(_RiGas);
  
  // needed for beta coefficient
  _mmasses.resize(nbSpecies);
  _library->getMolarMasses(_mmasses);
  
  _fcoeff.resize(nbSpecies);
  
  vector<CFuint> moleculeIDs;
  _library->setMoleculesIDs(moleculeIDs);
  vector<bool> flag(nbSpecies, false);
  for (CFuint i = 0; i < moleculeIDs.size(); ++i) {
    flag[moleculeIDs[i]] = true;
  }
  
  for (CFuint i = 0; i < nbSpecies; ++i) {
    _fcoeff[i] = (flag[i]) ? 2.5 : 1.5;
  }
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm::computePerturbedPhysicalData(const Framework::State& state,
						  const RealVector& pdataBkp,
						  RealVector& pdata,
						  CFuint iVar)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DNEQSymm::computePerturbedPhysicalData()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm::computeProjectedJacobian(const RealVector& normal, RealMatrix& jacob)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  cf_assert(getModel()->getNbScalarVars(1) == 1);
  
  const RealVector& lData = getModel()->getPhysicalData();
  SafePtr<PhysicalChemicalLibrary::ExtraData> eData = _library->getExtraData();
  
  const CFreal T = lData[EulerTerm::T];
  const CFreal rho = lData[EulerTerm::RHO];
  const CFreal p = lData[EulerTerm::P];
  const CFreal cvTr = eData->dEdT;
  const CFreal beta = p/(rho*T)/cvTr;
  
  cf_assert(_ys.size() == nbSpecies); 
  cf_assert(_alpha.size() == nbSpecies);
  
  if (p <= 0) {
    cout << "_ys = " << _ys << endl;
    cout << "rho = " << rho << endl;
    cout << "T   = " << T << endl;
    cout << "p   = " << p << endl;
    cout << "cvTr= " << cvTr << endl << endl;
    cf_assert(p > 0.0);
  }
  
  cf_assert(T > 0.0);
  cf_assert(rho > 0.0);
  
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal u = lData[EulerTerm::VX];
  const CFreal v = lData[EulerTerm::VY];
  const CFreal a2 = (1. + beta)*p/rho;
  const CFreal a = sqrt(a2);

  const CFreal U = u*nx + v*ny;
  ////////////////////////
  const CFuint uID  = nbSpecies;
  const CFuint vID  = nbSpecies+1;
  const CFuint eID  = nbSpecies+2;
  const CFuint evID = nbSpecies+3;
  ///////////////////////

  for (CFuint is = 0; is < nbSpecies; ++is) {
    for (CFuint js = 0; js < nbSpecies; ++js) {
      const CFreal Udelta = (js != is) ? 0.0 : U; 
      jacob(is,js) = Udelta;
    }
    jacob(is,uID) = 0.0;
    jacob(is,vID) = 0.0;
    jacob(is,eID) = 0.0;
    jacob(is,evID) = 0.0;
    
    jacob(uID,is) = 0.0;
    jacob(vID,is) = 0.0;
    jacob(eID,is) = 0.0;
    jacob(evID,is) = 0.0;
  }
  
  jacob(uID,uID) = U;
  jacob(uID,vID) = 0.0;
  jacob(uID,eID) = a*nx;
  jacob(uID,evID) = 0.0;
  
  jacob(vID,uID) = 0.0;
  jacob(vID,vID) = U;
  jacob(vID,eID) = a*ny;
  jacob(vID,evID) = 0.0;
  
  jacob(eID,uID) = a*nx;
  jacob(eID,vID) = a*ny;
  jacob(eID,eID) = U;
  jacob(eID,evID) = 0.0;
  
  jacob(evID,uID) = 0.0;
  jacob(evID,vID) = 0.0;
  jacob(evID,eID) = 0.0;
  jacob(evID,evID) = U;//jacob checked on 2009/4/23.
}


//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm::computeScalarJacobian(const RealVector& normal,
					   RealVector& jacob)
{
  const RealVector& lData = getModel()->getPhysicalData();
  cf_assert(jacob.size() == getModel()->getNbScalarVars(0) || jacob.size() == (getModel()->getNbScalarVars(0)+1));
  jacob = lData[EulerTerm::VX]*normal[XX] + lData[EulerTerm::VY]*normal[YY];
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm::splitJacobian_FullDecoupling(RealMatrix& jacobPlus,
						  RealMatrix& jacobMin,
						  RealVector& eValues,
						  const RealVector& normal)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  cf_assert(getModel()->getNbScalarVars(1) == 1);
  
  const RealVector& lData = getModel()->getPhysicalData();
  SafePtr<PhysicalChemicalLibrary::ExtraData> eData = _library->getExtraData();
    
  const CFreal T = lData[EulerTerm::T];
  const CFreal rho = lData[EulerTerm::RHO];
  const CFreal p = lData[EulerTerm::P];
  const CFreal cvTr = eData->dEdT;
  //const CFreal beta = p/(rho*T)/cvTr;
  
  const CFuint firstSpecies = getModel()->getFirstScalarVar(0);
  
  CFreal numBeta = 0.;
  CFreal denBeta = 0.;
  const CFuint start = (_library->presenceElectron()) ? 1 : 0;
  //Note: The definition of beta should NEVER include a contribution from the electrons (see in Gnoffo,Gupta,Shinn).
  
  for (CFuint i = start; i < nbSpecies; ++i) {
    const CFreal sigmai = lData[firstSpecies + i]/_mmasses[i];
    numBeta += sigmai;
    denBeta += sigmai*_fcoeff[i];
  }
  const CFreal beta = numBeta/denBeta;

  if (p <= 0) {
    cout << "_ys = " << _ys << endl;
    cout << "rho = " << rho << endl;
    cout << "T   = " << T << endl;
    cout << "p   = " << p << endl;
    cout << "cvTr= " << cvTr << endl << endl;
    cf_assert(p > 0.0);
  }
  
  cf_assert(T > 0.0);
  cf_assert(rho > 0.0);
        
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal u = lData[EulerTerm::VX];
  const CFreal v = lData[EulerTerm::VY];
  const CFreal V2 = lData[EulerTerm::V]*lData[EulerTerm::V];
  const CFreal q = 0.5*V2;
  const CFreal bq = beta*q;  
  const CFreal a2 = (1. + beta)*p/rho;
  const CFreal a = sqrt(a2);
  const CFreal U = u*nx + v*ny;
  
  cf_assert(_ys.size() == nbSpecies); 
  cf_assert(_alpha.size() == nbSpecies);
  
  for (CFuint is = 0; is < nbSpecies; ++is) {
    _ys[is] = lData[firstSpecies + is];
  }
  
  CFreal phi = 0.0; 
  if (_library->presenceElectron()) {
    // assume that electrons have 0 as mixture ID  
    phi = _RiGas[0]*_ys[0]/eData->dEvTv - beta;
  }
  else {
    phi =-beta;
  }
  
  CFreal Tq = 0.0;
  for (CFuint is = 0; is < nbSpecies; ++is) {
    // here T must be substituted with Tv if there is ionization
    // here phi has to be changed for the ionized case
    if (!_library->presenceElectron()) {
      Tq = T;
    }
    else {
      if (is == 0) {
	Tq = 0.0; //To implement: get vibrational temperature
      }
      else{
	Tq=T;
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
    //The previous commented portion of code shows differences up to 4e-10.
    //Thus we comment out the next assertion:
    //cf_assert(_RiGas[is]*Tq + bq - beta*(eData->energyTr)[is] == _alpha[is]);
  }
  
  const CFreal lambda1 = U;
  const CFreal lambda2 = U+a;
  const CFreal lambda3 = U-a;
  
  const CFreal lambda1P = max(0.0, lambda1);
  const CFreal lambda2P = max(0.0, lambda2);
  const CFreal lambda3P = max(0.0, lambda3);
      
  const CFreal lambda1M = min(0.0, lambda1);
  const CFreal lambda2M = min(0.0, lambda2);
  const CFreal lambda3M = min(0.0, lambda3);
  
  const CFreal lambda2P_Plus_Lambda3P = lambda2P+lambda3P;
  const CFreal lambda2P_Minus_Lambda3P = lambda2P-lambda3P;

  const CFreal lambda2M_Plus_Lambda3M = lambda2M+lambda3M;
  const CFreal lambda2M_Minus_Lambda3M = lambda2M-lambda3M;
  
  //  const CFreal lx = -ny; //Components of a vector normal to (nx , ny)
  // const CFreal ly = nx;
  
  //jacobPlus
  jacobPlus(0,0) = lambda1P*ny*ny+.5*lambda2P_Plus_Lambda3P*nx*nx;
//lambda1P*lx*lx+.5*lambda2P_Plus_Lambda3P*nx*nx;
  jacobPlus(0,1) = nx*ny*(.5*lambda2P_Plus_Lambda3P-lambda1P);
  jacobPlus(0,2) = .5*lambda2P_Minus_Lambda3P*nx;
  
  jacobPlus(1,0) = nx*ny*(.5*lambda2P_Plus_Lambda3P-lambda1P);
  jacobPlus(1,1) = lambda1P*nx*nx+.5*lambda2P_Plus_Lambda3P*ny*ny;
  jacobPlus(1,2) = .5*lambda2P_Minus_Lambda3P*ny;
  
  jacobPlus(2,0) = .5*lambda2P_Minus_Lambda3P*nx;
  jacobPlus(2,1) = .5*lambda2P_Minus_Lambda3P*ny;
  jacobPlus(2,2) = .5*lambda2P_Plus_Lambda3P;
  
  //jacobMin
  jacobMin(0,0) = lambda1M*ny*ny+.5*lambda2M_Plus_Lambda3M*nx*nx;
  jacobMin(0,1) = nx*ny*(.5*lambda2M_Plus_Lambda3M-lambda1M);
  jacobMin(0,2) = .5*lambda2M_Minus_Lambda3M*nx;
  
  jacobMin(1,0) = nx*ny*(.5*lambda2M_Plus_Lambda3M-lambda1M);
  jacobMin(1,1) = lambda1M*nx*nx+.5*lambda2M_Plus_Lambda3M*ny*ny;
  jacobMin(1,2) = .5*lambda2M_Minus_Lambda3M*ny;
  
  jacobMin(2,0) = .5*lambda2M_Minus_Lambda3M*nx;
  jacobMin(2,1) = .5*lambda2M_Minus_Lambda3M*ny;
  jacobMin(2,2) = .5*lambda2M_Plus_Lambda3M;
  
  eValues[0] = U;
  eValues[1] = U+a;
  eValues[2] = U-a;
  
  //debugging code. Everything worked OK.
  /*
  //tempPlus
  RealMatrix tempPlus(3,3);
  tempPlus(0,0) = lambda1P*lx*lx+.5*lambda2P_Plus_Lambda3P*nx*nx;
  tempPlus(0,1) = lambda1P*lx*ly+.5*lambda2P_Plus_Lambda3P*nx*ny;
  tempPlus(0,2) = .5*lambda2P_Minus_Lambda3P*nx;

  tempPlus(1,0) = lambda1P*lx*ly+.5*lambda2P_Plus_Lambda3P*nx*ny;
  tempPlus(1,1) = lambda1P*ly*ly+.5*lambda2P_Plus_Lambda3P*ny*ny;
  tempPlus(1,2) = .5*lambda2P_Minus_Lambda3P*ny;

  tempPlus(2,0) = .5*lambda2P_Minus_Lambda3P*nx;
  tempPlus(2,1) = .5*lambda2P_Minus_Lambda3P*ny;
  tempPlus(2,2) = .5*lambda2P_Plus_Lambda3P;

   //tempMin
   RealMatrix tempMin(3,3);
  tempMin(0,0) = lambda1M*lx*lx+.5*lambda2M_Plus_Lambda3M*nx*nx;
  tempMin(0,1) = lambda1M*lx*ly+.5*lambda2M_Plus_Lambda3M*nx*ny;
  tempMin(0,2) = .5*lambda2M_Minus_Lambda3M*nx;

  tempMin(1,0) = lambda1M*lx*ly+.5*lambda2M_Plus_Lambda3M*nx*ny;
  tempMin(1,1) = lambda1M*ly*ly+.5*lambda2M_Plus_Lambda3M*ny*ny;
  tempMin(1,2) = .5*lambda2M_Minus_Lambda3M*ny;

  tempMin(2,0) = .5*lambda2M_Minus_Lambda3M*nx;
  tempMin(2,1) = .5*lambda2M_Minus_Lambda3M*ny;
  tempMin(2,2) = .5*lambda2M_Plus_Lambda3M;
  
    RealMatrix mat(3,3);
    mat= tempPlus+tempMin;
    
    cout << endl;
    cout <<"tempPlus" << endl;
    cout << tempPlus << endl << endl;
    
    cout << endl;
    cout <<"tempMin" << endl;
    cout << tempMin << endl << endl;

    cout << endl;
    cout <<"tempPlus+tempMin" << endl;
    cout << mat << endl << endl;

    
    RealMatrix mat2(3,3);
    mat2(0,0) = U;
    mat2(0,1) = 0.0;
    mat2(0,2) = a*nx;
    
    mat2(1,0) = 0.0;
    mat2(1,1) = U;
    mat2(1,2) = a*ny;
    
    mat2(2,0) = a*nx;
    mat2(2,1) = a*ny;
    mat2(2,2) = U;
    cout  << endl;
    cout << "jacob"<< endl << mat2 << endl << endl;
    
    mat =mat-mat2;
    cout  << endl;
    cout << "tempPlus+tempMin-jacob"<< endl <<mat << endl << endl;
    */

  //degugging infos
  CFLogDebugMax( "RightEigenvectors @Euler2DNEQSymm::splitJacobian_FullDecoupling" << "\n"
		 << _rightEv << "\n");
  CFLogDebugMax( "LeftEigenvectors @Euler2DNEQSymm::splitJacobian_FullDecoupling" << "\n"
		 << _leftEv << "\n");
  CFLogDebugMax( "EigenValues @Euler2DNEQSymm::splitJacobian_FullDecoupling" << "\n"
		 << eValues << "\n" << "\n"); 
  //
  //  cout<<"Aborting from method splitJacobian in the Object Euler2DNEQSymm."<<endl;
  // abort();
  //Comment back, from here
  /*for (CFuint count=0;count<lData.size();count++){
    cout << endl<< count << " index => " << lData[count] << endl;
  }
  cout<<endl<<"Aborted from method splitJacobian"<<endl;
  abort();*/
  //to here
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm::splitJacobian_PartialDecoupling(RealMatrix& jacobPlus,
						     RealMatrix& jacobMin,
						     RealVector& eValues,
						     const RealVector& normal)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  cf_assert(getModel()->getNbScalarVars(1) == 1);
  
  const RealVector& lData = getModel()->getPhysicalData();
  SafePtr<PhysicalChemicalLibrary::ExtraData> eData = _library->getExtraData();
    
  const CFreal T = lData[EulerTerm::T];
  const CFreal rho = lData[EulerTerm::RHO];
  const CFreal p = lData[EulerTerm::P];
  const CFreal cvTr = eData->dEdT;
  //const CFreal beta = p/(rho*T)/cvTr;
  
  const CFuint firstSpecies = getModel()->getFirstScalarVar(0);

  CFreal numBeta = 0.;
  CFreal denBeta = 0.;
  const CFuint start = (_library->presenceElectron()) ? 1 : 0;
  //Note: The definition of beta should NEVER include a contribution from the electrons (see in Gnoffo,Gupta,Shinn).
  
  for (CFuint i = start; i < nbSpecies; ++i) {
    const CFreal sigmai = lData[firstSpecies + i]/_mmasses[i];
    numBeta += sigmai;
    denBeta += sigmai*_fcoeff[i];
  }
  const CFreal beta = numBeta/denBeta;
       
  if (p <= 0) {
    cout << "_ys = " << _ys << endl;
    cout << "rho = " << rho << endl;
    cout << "T   = " << T << endl;
    cout << "p   = " << p << endl;
    cout << "cvTr= " << cvTr << endl << endl;
    cf_assert(p > 0.0);
  }
  
  cf_assert(T > 0.0);
  cf_assert(rho > 0.0);
        
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal u = lData[EulerTerm::VX];
  const CFreal v = lData[EulerTerm::VY];
  const CFreal V2 = lData[EulerTerm::V]*lData[EulerTerm::V];
  const CFreal q = 0.5*V2;
  const CFreal bq = beta*q;  
  const CFreal a2 = (1. + beta)*p/rho;
  const CFreal a = sqrt(a2);
  //const CFreal a = lData[EulerTerm::A];
  //const CFreal a2 = lData[EulerTerm::A]*lData[EulerTerm::A];
  
  const CFreal U = u*nx + v*ny;
    
  cf_assert(_ys.size() == nbSpecies); 
  cf_assert(_alpha.size() == nbSpecies);
  
  for (CFuint is = 0; is < nbSpecies; ++is) {
    _ys[is] = lData[firstSpecies + is];
  }
  
  CFreal phi = 0.0; 
  if (_library->presenceElectron()) {
    // assume that electrons have 0 as mixture ID  
    phi = _RiGas[0]*_ys[0]/eData->dEvTv - beta;
  }
  else {
    phi =-beta;
  }
  
  CFreal Tq = 0.0;
  for (CFuint is = 0; is < nbSpecies; ++is) {
    // here T must be substituted with Tv if there is ionization
    // here phi has to be changed for the ionized case
    if (!_library->presenceElectron()) {
      Tq = T;
    }
    else {
      if (is == 0) {
	Tq = 0.0; //To implement: get vibrational temperature
      }
      else{
	Tq=T;
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
    //The previous commented portion of code shows differences up to 4e-10.
    //Thus we comment out the next assertion:
    //cf_assert(_RiGas[is]*Tq + bq - beta*(eData->energyTr)[is] == _alpha[is]);
  }
  
  //_rightEv = 0.0;
  //_leftEv = 0.0;
  const CFreal lambda1 = U;
  const CFreal lambda2 = U+a;
  const CFreal lambda3 = U-a;
  
  const CFreal lambda1P = max(0.0, lambda1);
  const CFreal lambda2P = max(0.0, lambda2);
  const CFreal lambda3P = max(0.0, lambda3);
      
  const CFreal lambda1M = min(0.0, lambda1);
  const CFreal lambda2M = min(0.0, lambda2);
  const CFreal lambda3M = min(0.0, lambda3);
  
  const CFreal lambda2P_Plus_Lambda3P = lambda2P+lambda3P;
  const CFreal lambda2P_Minus_Lambda3P = lambda2P-lambda3P;

  const CFreal lambda2M_Plus_Lambda3M = lambda2M+lambda3M;
  const CFreal lambda2M_Minus_Lambda3M = lambda2M-lambda3M;
  
  const CFreal lx = -ny; //Components of a vector normal to (nx , ny)
  const CFreal ly = nx;
  
  //jacobPlus
  
  jacobPlus(0,0) = lambda1P*lx*lx+.5*lambda2P_Plus_Lambda3P*nx*nx;
  jacobPlus(0,1) = lambda1P*lx*ly+.5*lambda2P_Plus_Lambda3P*nx*ny;
  jacobPlus(0,2) = .5*lambda2P_Minus_Lambda3P*nx;
  jacobPlus(0,3) = 0.0;
  
  jacobPlus(1,0) = lambda1P*lx*ly+.5*lambda2P_Plus_Lambda3P*nx*ny;
  jacobPlus(1,1) = lambda1P*ly*ly+.5*lambda2P_Plus_Lambda3P*ny*ny;
  jacobPlus(1,2) = .5*lambda2P_Minus_Lambda3P*ny;
  jacobPlus(1,3) = 0.0;
  
  jacobPlus(2,0) = .5*lambda2P_Minus_Lambda3P*nx;
  jacobPlus(2,1) = .5*lambda2P_Minus_Lambda3P*ny;
  jacobPlus(2,2) = .5*lambda2P_Plus_Lambda3P;
  jacobPlus(2,3) = 0.0;
  
  jacobPlus(3,0) = 0.0;
  jacobPlus(3,1) = 0.0;
  jacobPlus(3,2) = 0.0;
  jacobPlus(3,3) = lambda1P;
  
  //jacobMin
  
  jacobMin(0,0) = lambda1M*lx*lx+.5*lambda2M_Plus_Lambda3M*nx*nx;
  jacobMin(0,1) = lambda1M*lx*ly+.5*lambda2M_Plus_Lambda3M*nx*ny;
  jacobMin(0,2) = .5*lambda2M_Minus_Lambda3M*nx;
  jacobMin(0,3) = 0.0;
  
  jacobMin(1,0) = lambda1M*lx*ly+.5*lambda2M_Plus_Lambda3M*nx*ny;
  jacobMin(1,1) = lambda1M*ly*ly+.5*lambda2M_Plus_Lambda3M*ny*ny;
  jacobMin(1,2) = .5*lambda2M_Minus_Lambda3M*ny;
  jacobMin(1,3) = 0.0;
  
  jacobMin(2,0) = .5*lambda2M_Minus_Lambda3M*nx;
  jacobMin(2,1) = .5*lambda2M_Minus_Lambda3M*ny;
  jacobMin(2,2) = .5*lambda2M_Plus_Lambda3M;
  jacobMin(2,3) = 0.0;
  
  jacobMin(3,0) = 0.0;
  jacobMin(3,1) = 0.0;
  jacobMin(3,2) = 0.0;
  jacobMin(3,3) = lambda1M;
  
  eValues[0] = U;
  eValues[1] = U+a;
  eValues[2] = U-a;
  eValues[3] = U;
  
  //debugging code. Everything worked OK.
  /*RealMatrix mat(jacobSize, jacobSize);
    mat= jacobPlus+jacobMin;
    
    cout << endl;
    cout <<"jacobPlus" << endl;
    cout << jacobPlus << endl << endl;
    
    cout << endl;
    cout <<"jacobMin" << endl;
    cout << jacobMin << endl << endl;

    cout << endl;
    cout <<"jacobPlus+jacobMin" << endl;
    cout << mat << endl << endl;

    
    RealMatrix mat2(jacobSize, jacobSize);
    mat2(0,0) = U;
    mat2(0,1) = 0.0;
    mat2(0,2) = a*nx;
    mat2(0,3) = 0.0;

    mat2(1,0) = 0.0;
    mat2(1,1) = U;
    mat2(1,2) = a*ny;
    mat2(1,3) = 0.0;

    mat2(2,0) = a*nx;
    mat2(2,1) = a*ny;
    mat2(2,2) = U;
    mat2(2,3) = 0.0;
    
    mat2(3,0) = 0.0;
    mat2(3,1) = 0.0;
    mat2(3,2) = 0.0;
    mat2(3,3) = U;
    
    mat =mat-mat2;
    cout  << endl;
    cout << "jacobPlus+jacobMin-jacob"<< endl <<mat << endl << endl;
    */

  //degugging infos
  CFLogDebugMax( "RightEigenvectors @Euler2DNEQSymm::splitJacobian_PartialDecoupling" << "\n"
		 << _rightEv << "\n");
  CFLogDebugMax( "LeftEigenvectors @Euler2DNEQSymm::splitJacobian_PartialDecoupling" << "\n"
		 << _leftEv << "\n");
  CFLogDebugMax( "EigenValues @Euler2DNEQSymm::splitJacobian_PartialDecoupling" << "\n"
		 << eValues << "\n" << "\n"); 
  //
  //  cout<<"Aborting from method splitJacobian in the Object Euler2DNEQSymm."<<endl;
  // abort();
  //Comment back, from here
  /*for (CFuint count=0;count<lData.size();count++){
    cout << endl<< count << " index => " << lData[count] << endl;
  }
  cout<<endl<<"Aborted from method splitJacobian"<<endl;
  abort();*/
  //to here
  
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm::splitJacobian_NoDecoupling(RealMatrix& jacobPlus,
						RealMatrix& jacobMin,
						RealVector& eValues,
						const RealVector& normal)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  cf_assert(getModel()->getNbScalarVars(1) == 1);
  
  const RealVector& lData = getModel()->getPhysicalData();
  SafePtr<PhysicalChemicalLibrary::ExtraData> eData = _library->getExtraData();
    
  const CFreal T = lData[EulerTerm::T];
  const CFreal rho = lData[EulerTerm::RHO];
  const CFreal p = lData[EulerTerm::P];
  const CFreal cvTr = eData->dEdT;
  //const CFreal beta = p/(rho*T)/cvTr;
  
  const CFuint firstSpecies = getModel()->getFirstScalarVar(0);

  CFreal numBeta = 0.;
  CFreal denBeta = 0.;
  const CFuint start = (_library->presenceElectron()) ? 1 : 0;
  //Note: The definition of beta should NEVER include a contribution from the electrons (see in Gnoffo,Gupta,Shinn).
  
  for (CFuint i = start; i < nbSpecies; ++i) {
    const CFreal sigmai = lData[firstSpecies + i]/_mmasses[i];
    numBeta += sigmai;
    denBeta += sigmai*_fcoeff[i];
  }
  const CFreal beta = numBeta/denBeta;
  
  if (p <= 0) {
    cout << "_ys = " << _ys << endl;
    cout << "rho = " << rho << endl;
    cout << "T   = " << T << endl;
    cout << "p   = " << p << endl;
    cout << "cvTr= " << cvTr << endl << endl;
    cf_assert(p > 0.0);
  }
  
  cf_assert(T > 0.0);
  cf_assert(rho > 0.0);
        
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal u = lData[EulerTerm::VX];
  const CFreal v = lData[EulerTerm::VY];
  const CFreal V2 = lData[EulerTerm::V]*lData[EulerTerm::V];
  const CFreal q = 0.5*V2;
  const CFreal bq = beta*q;  
  const CFreal a2 = (1. + beta)*p/rho;
  const CFreal a = sqrt(a2);
  const CFreal ova = 1./a;
  const CFreal U = u*nx + v*ny;
  
  const CFuint nbSp = nbSpecies;
  const CFuint nbSpPlus1 = nbSpecies+1; 
  const CFuint nbSpPlus2 = nbSpecies+2; 
  const CFuint nbSpPlus3 = nbSpecies+3;
  
  const CFuint uID  = nbSpecies;
  const CFuint vID  = nbSpecies+1;
  const CFuint eID  = nbSpecies+2;
  const CFuint evID = nbSpecies+3;
  
  cf_assert(_ys.size() == nbSpecies); 
  cf_assert(_alpha.size() == nbSpecies);
  
  for (CFuint is = 0; is < nbSpecies; ++is) {
    _ys[is] = lData[firstSpecies + is];
  }
  
  CFreal phi = 0.0; 
  if (_library->presenceElectron()) {
    // assume that electrons have 0 as mixture ID  
    phi = _RiGas[0]*_ys[0]/eData->dEvTv - beta;
  }
  else {
    phi =-beta;
  }
  
  CFreal Tq = 0.0;
  for (CFuint is = 0; is < nbSpecies; ++is) {
    // here T must be substituted with Tv if there is ionization
    // here phi has to be changed for the ionized case
    if (!_library->presenceElectron()) {
      Tq = T;
    }
    else {
      if (is == 0) {
	Tq = 0.0; //To implement: get vibrational temperature
      }
      else{
	Tq=T;
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
    //The previous commented portion of code shows differences up to 4e-10.
    //Thus we comment out the next assertion:
    //cf_assert(_RiGas[is]*Tq + bq - beta*(eData->energyTr)[is] == _alpha[is]);
  }
  
  for (CFuint is = 0; is < nbSpecies; ++is) {
    for (CFuint js = 0; js < nbSpecies; ++js) {
      const CFreal minusRhoDelta = (js != is) ? 0.0 : -rho; 
      _rightEv(is, js) = minusRhoDelta;
    }
    
    _rightEv(is,nbSp) = 0.0;
    _rightEv(is,nbSpPlus1) = 0.0;
    _rightEv(is,nbSpPlus2) = 0.0;
    _rightEv(is,nbSpPlus3) = 0.0;
    
    
    _rightEv(nbSp,is) = 0.0;
    _rightEv(nbSpPlus1,is) = 0.0;
    _rightEv(nbSpPlus2,is) = 0.0;
    _rightEv(nbSpPlus3,is) = 0.0;
  }  
  
  _rightEv(nbSp,nbSp) = -ny;
  _rightEv(nbSp,nbSpPlus1) = 0.5*ova*nx;
  _rightEv(nbSp,nbSpPlus2) = -0.5*ova*nx;
  _rightEv(nbSp,nbSpPlus3) = 0.0;
  
  _rightEv(nbSpPlus1,nbSp) = nx;
  _rightEv(nbSpPlus1,nbSpPlus1) = 0.5*ova*ny;
  _rightEv(nbSpPlus1,nbSpPlus2) = -0.5*ova*ny;
  _rightEv(nbSpPlus1,nbSpPlus3) = 0.0;  
  
  _rightEv(nbSpPlus2,nbSp) = 0.0;
  _rightEv(nbSpPlus2,nbSpPlus1) = 0.5*ova;
  _rightEv(nbSpPlus2,nbSpPlus2) = 0.5*ova;
  _rightEv(nbSpPlus2,nbSpPlus3) = 0.0;
  
  _rightEv(nbSpPlus3,nbSp) = 0.0;
  _rightEv(nbSpPlus3,nbSpPlus1) = 0.0;
  _rightEv(nbSpPlus3,nbSpPlus2) = 0.0;
  _rightEv(nbSpPlus3,nbSpPlus3) = -1.0;
  
  _rightEv /= rho;//_rightEv checked on 2009/3/30.
  //
  //cout << "R = "<< endl;
  //cout << _rightEv << endl <<endl;
  //
  
  for (CFuint is = 0; is < nbSpecies; ++is) {
    for (CFuint js = 0; js < nbSpecies; ++js) {
      const CFreal minusdelta = (js != is) ? 0.0 : -1.0; 
      _leftEv(is,js) = minusdelta;
    }
    
    _leftEv(is,nbSp) = 0.0;
    _leftEv(is,nbSpPlus1) = 0.0;
    _leftEv(is,nbSpPlus2) = 0.0;
    _leftEv(is,nbSpPlus3) = 0.0;
    
    _leftEv(nbSp,is) = 0.0;
    _leftEv(nbSpPlus1,is) = 0.0;
    _leftEv(nbSpPlus2,is) = 0.0;
    _leftEv(nbSpPlus3,is) = 0.0;
  }  
  
  _leftEv(nbSp,nbSp) = -rho*ny;
  _leftEv(nbSp,nbSpPlus1) = rho*nx;
  _leftEv(nbSp,nbSpPlus2) = 0.0;
  _leftEv(nbSp,nbSpPlus3) = 0.0;
  
  
  _leftEv(nbSpPlus1,nbSp) = rho*a*nx;
  _leftEv(nbSpPlus1,nbSpPlus1) = rho*a*ny;
  _leftEv(nbSpPlus1,nbSpPlus2) = rho*a;
  _leftEv(nbSpPlus1,nbSpPlus3) = 0.0;
  
  _leftEv(nbSpPlus2,nbSp) = -rho*a*nx;
  _leftEv(nbSpPlus2,nbSpPlus1) = -rho*a*ny;
  _leftEv(nbSpPlus2,nbSpPlus2) = rho*a;
  _leftEv(nbSpPlus2,nbSpPlus3) = 0.0;
  
  _leftEv(nbSpPlus3,nbSp) = 0.0; 
  _leftEv(nbSpPlus3,nbSpPlus1) = 0.0;
  _leftEv(nbSpPlus3,nbSpPlus2) = 0.0;
  _leftEv(nbSpPlus3,nbSpPlus3) = -rho;//_leftEv checked on 2009/3/30.
  
  //Comment back from here to
  //cout << "L = "<< endl;
  //cout << _leftEv << endl <<endl;
  
  //RealMatrix mat(_rightEv.nbRows(), _leftEv.nbRows());
  //mat= _rightEv*_leftEv;
  //cout <<"R*L" << endl;
  //cout << mat << endl << endl;
  
  //mat = _leftEv*_rightEv;
  //cout <<"L*R" << endl;
  //cout << mat << endl << endl;// matrices _leftEv & _rightEv are the inverse of the other.
  //here
  
  //Beginning of debugging code.
  
  //End of debugging code.
  
  const CFreal lambda1 = U;
  const CFreal lambda2 = U+a;
  const CFreal lambda3 = U-a;
  
  const CFreal lambda1P = max(0.0, lambda1);
  const CFreal lambda2P = max(0.0, lambda2);
  const CFreal lambda3P = max(0.0, lambda3);
      
  const CFreal lambda1M = min(0.0, lambda1);
  const CFreal lambda2M = min(0.0, lambda2);
  const CFreal lambda3M = min(0.0, lambda3);
  
  const CFreal lambda2P_Plus_Lambda3P = lambda2P+lambda3P;
  const CFreal lambda2P_Minus_Lambda3P = lambda2P-lambda3P;

  const CFreal lambda2M_Plus_Lambda3M = lambda2M+lambda3M;
  const CFreal lambda2M_Minus_Lambda3M = lambda2M-lambda3M;
      
  eValues = U;
  eValues[nbSpPlus1] += a;
  eValues[nbSpPlus2] -= a;

  // compute the eigen values + and -
  if (_jacobDissip > 0.0) {
    // modified eigenvalues to cure carbuncle
    const CFreal j2 = _jacobDissip*_jacobDissip;
    for (CFuint iEq = 0; iEq < eValues.size(); ++iEq) {
      _eValuesP[iEq] = max(0.,eValues[iEq]);
      const CFreal evP = _eValuesP[iEq];
      _eValuesP[iEq] = 0.5*(evP + sqrt(evP*evP + j2*a2));
      
      _eValuesM[iEq] = min(0.,eValues[iEq]);
      const CFreal evM = _eValuesM[iEq];
      _eValuesM[iEq] = 0.5*(evM - sqrt(evM*evM + j2*a2));
    }
  }
  else {
    for (CFuint iEq = 0; iEq < eValues.size(); ++iEq) {
      _eValuesP[iEq] = max(0.,eValues[iEq]);
      _eValuesM[iEq] = min(0.,eValues[iEq]);
    }
  }
  
  // compute jacobian + and -
  
  //jacobPlus = _rightEv*(_eValuesP*_leftEv);
  for (CFuint is = 0; is < nbSpecies; ++is) {
    for (CFuint js = 0; js < nbSpecies; ++js) {
      jacobPlus(is,js) =  (js != is) ? 0.0 : lambda1P; 
    }
    jacobPlus(is,uID) = 0.0;
    jacobPlus(is,vID) = 0.0;
    jacobPlus(is,eID) = 0.0;
    jacobPlus(is,evID) = 0.0;
    
    jacobPlus(uID,is) = 0.0;
    jacobPlus(vID,is) = 0.0;
    jacobPlus(eID,is) = 0.0;
    jacobPlus(evID,is) = 0.0;
  }
  
  jacobPlus(uID,uID) = 0.5*lambda2P_Plus_Lambda3P*nx*nx+lambda1P*ny*ny;
  jacobPlus(uID,vID) = nx*ny*(0.5*lambda2P_Plus_Lambda3P-lambda1P);
  jacobPlus(uID,eID) = 0.5*lambda2P_Minus_Lambda3P*nx;
  jacobPlus(uID,evID) = 0.0;
  
  jacobPlus(vID,uID) = nx*ny*(0.5*lambda2P_Plus_Lambda3P-lambda1P);
  jacobPlus(vID,vID) = lambda1P*nx*nx+0.5*lambda2P_Plus_Lambda3P*ny*ny;
  jacobPlus(vID,eID) = 0.5*lambda2P_Minus_Lambda3P*ny;
  jacobPlus(vID,evID) = 0.0;
  
  jacobPlus(eID,uID) = 0.5*lambda2P_Minus_Lambda3P*nx;
  jacobPlus(eID,vID) = 0.5*lambda2P_Minus_Lambda3P*ny;
  jacobPlus(eID,eID) = 0.5*lambda2P_Plus_Lambda3P;
  jacobPlus(eID,evID) = 0.0;
  
  jacobPlus(evID,uID) = 0.0;
  jacobPlus(evID,vID) = 0.0;
  jacobPlus(evID,eID) = 0.0;
  jacobPlus(evID,evID) = lambda1P;
  
  //jacobMin  = _rightEv*(_eValuesM*_leftEv);      
  for (CFuint is = 0; is < nbSpecies; ++is) {
    for (CFuint js = 0; js < nbSpecies; ++js) {
      jacobMin(is,js) =  (js != is) ? 0.0 : lambda1M; 
    }
    jacobMin(is,uID) = 0.0;
    jacobMin(is,vID) = 0.0;
    jacobMin(is,eID) = 0.0;
    jacobMin(is,evID) = 0.0;
    
    jacobMin(uID,is) = 0.0;
    jacobMin(vID,is) = 0.0;
    jacobMin(eID,is) = 0.0;
    jacobMin(evID,is) = 0.0;
  }
  
  jacobMin(uID,uID) = 0.5*lambda2M_Plus_Lambda3M*nx*nx+lambda1M*ny*ny;
  jacobMin(uID,vID) = nx*ny*(0.5*lambda2M_Plus_Lambda3M-lambda1M);
  jacobMin(uID,eID) = 0.5*lambda2M_Minus_Lambda3M*nx;
  jacobMin(uID,evID) = 0.0;
  
  jacobMin(vID,uID) = nx*ny*(0.5*lambda2M_Plus_Lambda3M-lambda1M);
  jacobMin(vID,vID) = lambda1M*nx*nx+0.5*lambda2M_Plus_Lambda3M*ny*ny;
  jacobMin(vID,eID) = 0.5*lambda2M_Minus_Lambda3M*ny;
  jacobMin(vID,evID) = 0.0;
  
  jacobMin(eID,uID) = 0.5*lambda2M_Minus_Lambda3M*nx;
  jacobMin(eID,vID) = 0.5*lambda2M_Minus_Lambda3M*ny;
  jacobMin(eID,eID) = 0.5*lambda2M_Plus_Lambda3M;
  jacobMin(eID,evID) = 0.0;
  
  jacobMin(evID,uID) = 0.0;
  jacobMin(evID,vID) = 0.0;
  jacobMin(evID,eID) = 0.0;
  jacobMin(evID,evID) = lambda1M;
  
  
  //Debugging code for jacobPlus and jacobMin. Everything went OK.
  /*
    RealMatrix mat(_rightEv.nbRows(), _leftEv.nbRows());
    mat= jacobPlus+jacobMin;
    cout << endl;
    cout <<"jacobPlus+jacobMin" << endl;
  cout << mat << endl << endl;
  
  RealMatrix jacob(_rightEv.nbRows(), _leftEv.nbRows());
  computeProjectedJacobian( normal,  jacob);
  
  mat =mat-jacob;
  cout  << endl;
  cout << "jacobPlus+jacobMin-jacob"<< endl <<mat << endl << endl;
  */


  //degugging infos
  CFLogDebugMax( "RightEigenvectors @Euler2DNEQSymm::splitJacobian" << "\n"
		 << _rightEv << "\n");
  CFLogDebugMax( "LeftEigenvectors @Euler2DNEQSymm::splitJacobian" << "\n"
		 << _leftEv << "\n");
  CFLogDebugMax( "EigenValues @Euler2DNEQSymm::splitJacobian" << "\n"
		 << eValues << "\n" << "\n"); 
  //
  //  cout<<"Aborting from method splitJacobian in the Object Euler2DNEQSymm."<<endl;
  // abort();
  //Comment back, from here
  /*for (CFuint count=0;count<lData.size();count++){
    cout << endl<< count << " index => " << lData[count] << endl;
  }
  cout<<endl<<"Aborted from method splitJacobian"<<endl;
  abort();*/
  //to here  
}

/////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm::setStateVelocityIDs (std::vector<CFuint>& velIDs) 
{
  const CFuint nbSpecies = _library->getNbSpecies();
  velIDs.resize(2); velIDs[XX] = nbSpecies; velIDs[YY] = nbSpecies + 1; 
}

//////////////////////////////////////////////////////////////////////////////

    } // end of namespace NEQ

  } // end of namespace Physics

} // end of namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
