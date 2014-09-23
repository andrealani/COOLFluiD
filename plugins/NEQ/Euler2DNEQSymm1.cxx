#include "NEQ.hh"
#include "Euler2DNEQSymm1.hh"
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

Environment::ObjectProvider<Euler2DNEQSymm1, ConvectiveVarSet, NEQModule, 1>
euler2DNEQSymm1Provider("Euler2DNEQSymm1");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQSymm1::Euler2DNEQSymm1(Common::SafePtr<BaseTerm> term) :
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

Euler2DNEQSymm1::~Euler2DNEQSymm1()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm1::computeEigenValuesVectors(RealMatrix& rightEv,
                                           RealMatrix& leftEv,
                                           RealVector& eValues,
                                           const RealVector& normal)
{
  throw Common::NotImplementedException(FromHere(), "Euler2DNEQSymm1::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler2DNEQSymm1::getBlockSeparator() const
{
  return Framework::PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm1::splitJacobian(RealMatrix& jacobPlus,
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
    
  CFreal sumAlphaYs = 0.0;
     
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
  //const CFreal V = v*nx -u*ny;
  //const CFreal H = lData[EulerTerm::H];
  const CFuint firstTv = getModel()->getFirstScalarVar(1);
  const CFreal ev = lData[firstTv];
  const CFreal a2 = (1. + beta)*p/rho;
  const CFreal a = sqrt(a2);
  //const CFreal a = lData[EulerTerm::A];
  //const CFreal a2 = lData[EulerTerm::A]*lData[EulerTerm::A];
  const CFreal ova = 1./a;
  //const CFreal ova2 = 1./a2;

  const CFreal U = u*nx + v*ny;
  
  const CFuint nbSp = nbSpecies;
  const CFuint nbSpPlus1 = nbSpecies+1; 
  const CFuint nbSpPlus2 = nbSpecies+2; 
  const CFuint nbSpPlus3 = nbSpecies+3;
  ////////////////////////
  //const CFuint uID  = nbSpecies;
  //const CFuint vID  = nbSpecies+1;
  //const CFuint eID  = nbSpecies+2;
  //const CFuint evID = nbSpecies+3;
  /////////////////////////
  
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
    sumAlphaYs += _alpha[is]*_ys[is];   
  }
 
  //_rightEv = 0.0;
  //_leftEv = 0.0;

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
    _rightEv(nbSpPlus3,is) = ev;
  }  
  
  _rightEv(nbSp,nbSp) = -ny;
  _rightEv(nbSp,nbSpPlus1) = 0.5*nx*ova;
  _rightEv(nbSp,nbSpPlus2) = -0.5*nx*ova;
  _rightEv(nbSp,nbSpPlus3) = 0.0;

  _rightEv(nbSpPlus1,nbSp) = nx;
  _rightEv(nbSpPlus1,nbSpPlus1) = 0.5*ny*ova;
  _rightEv(nbSpPlus1,nbSpPlus2) = -0.5*ny*ova;
  _rightEv(nbSpPlus1,nbSpPlus3) = 0.0;  

  _rightEv(nbSpPlus2,nbSp) = 0.0;
  _rightEv(nbSpPlus2,nbSpPlus1) = 0.5*ova;
  _rightEv(nbSpPlus2,nbSpPlus2) = 0.5*ova;
  _rightEv(nbSpPlus2,nbSpPlus3) = 0.0;
  
  _rightEv(nbSpPlus3,nbSp) = 0.0;
  _rightEv(nbSpPlus3,nbSpPlus1) = 0.5*ev;
  _rightEv(nbSpPlus3,nbSpPlus2) = 0.5*ev;
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
    _leftEv(nbSpPlus3,is) = -ev;
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
  _leftEv(nbSpPlus3,nbSpPlus2) = rho*ev*a;
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
  jacobPlus = _rightEv*(_eValuesP*_leftEv);
  jacobMin  = _rightEv*(_eValuesM*_leftEv);
  

  //degugging infos
  CFLogDebugMax( "RightEigenvectors @Euler2DNEQSymm1::splitJacobian" << "\n"
		 << _rightEv << "\n");
  CFLogDebugMax( "LeftEigenvectors @Euler2DNEQSymm1::splitJacobian" << "\n"
		 << _leftEv << "\n");
  CFLogDebugMax( "EigenValues @Euler2DNEQSymm1::splitJacobian" << "\n"
		 << eValues << "\n" << "\n"); 
  //
  //  cout<<"Aborting from method splitJacobian in the Object Euler2DNEQSymm1."<<endl;
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

void Euler2DNEQSymm1::computePhysicalData(const State& state,
				     RealVector& data)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DNEQSymm1::computePhysicalData()");
}
      
//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm1::computeStateFromPhysicalData(const RealVector& data,
					 State& state)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DNEQSymm1::computeStateFromPhysicalData()");
  
  //  state[0] = data[EulerTerm::RHO];
  //   state[1] = data[EulerTerm::RHO]*data[EulerTerm::VX];
  //   state[2] = data[EulerTerm::RHO]*data[EulerTerm::VY];
  //   state[3] = data[EulerTerm::RHO]*data[EulerTerm::E];
  
  //   // Set the species
  //   const CFuint firstSpecies = getModel()->getFirstScalarVar(0);
  //   const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  
  //   for (CFuint ie = 0; ie < nbSpecies; ++ie){
  //     state[4 + ie] = data[EulerTerm::RHO]*data[firstSpecies + ie];
  //   }
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler2DNEQSymm1::getSpeed(const State& state) const
{
  const CFreal u = state[1]/state[0];
  const CFreal v = state[2]/state[0];
  return sqrt(u*u + v*v);
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm1::setDimensionalValues(const State& state,
                                          RealVector& result)
{
  throw Common::NotImplementedException
      (FromHere(), "Euler2DNEQSymm1::setDimensionalValues()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm1::setAdimensionalValues(const Framework::State& state,
                             RealVector& result)
{
  throw Common::NotImplementedException
      (FromHere(), "Euler2DNEQSymm1::setAdimensionalValues() not implemented");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm1::setDimensionalValuesPlusExtraValues
(const State& state, RealVector& result,
 RealVector& extra)
{
  throw Common::NotImplementedException
    (FromHere(), "Euler2DNEQSymm1::setDimensionalValuesPlusExtraValues()");
}
      
//////////////////////////////////////////////////////////////////////////////

vector<std::string> Euler2DNEQSymm1::getExtraVarNames() const
{
  cf_assert(_library.isNotNull());

  vector<std::string> names(3);
  names[0] = "rho";
  names[1] = "H";
  names[2] = "M";

  return names;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm1::setup()
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

void Euler2DNEQSymm1::computePerturbedPhysicalData(const Framework::State& state,
						  const RealVector& pdataBkp,
						  RealVector& pdata,
						  CFuint iVar)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DNEQSymm1::computePerturbedPhysicalData()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm1::computeProjectedJacobian(const RealVector& normal, RealMatrix& jacob)
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
  
  // const CFuint firstSpecies = getModel()->getFirstScalarVar(0);
  
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
  // const CFreal V2 = lData[EulerTerm::V]*lData[EulerTerm::V];
  // const CFreal q = 0.5*V2;
  // const CFreal V = v*nx -u*ny;
  // const CFreal H = lData[EulerTerm::H];
  const CFuint firstTv = getModel()->getFirstScalarVar(1);
  const CFreal ev = lData[firstTv];
  const CFreal a2 = (1. + beta)*p/rho;
  const CFreal a = sqrt(a2);

  const CFreal U = u*nx + v*ny;
  // const CFuint nbSpPlus1 = nbSpecies+1; 
  // const CFuint nbSpPlus2 = nbSpecies+2; 
  // const CFuint nbSpPlus3 = nbSpecies+3;
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
  
  jacob(evID,uID) = ev*a2*nx;
  jacob(evID,vID) = ev*a2*ny;
  jacob(evID,eID) = 0.0;
  jacob(evID,evID) = U;//jacob checked on 2009/4/23.
}


//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm1::setStateVelocityIDs (std::vector<CFuint>& velIDs) 
{
  const CFuint nbSpecies = _library->getNbSpecies();
  velIDs.resize(2); velIDs[XX] = nbSpecies; velIDs[YY] = nbSpecies + 1; 
}

//////////////////////////////////////////////////////////////////////////////

    } // end of namespace NEQ

  } // end of namespace Physics

} // end of namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
