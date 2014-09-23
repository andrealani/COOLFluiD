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
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint nbTv = getModel()->getNbScalarVars(1);
  cf_assert(nbTv == 1);

  const RealVector& lData = getModel()->getPhysicalData();
  SafePtr<PhysicalChemicalLibrary::ExtraData> eData = _library->getExtraData();

  const CFreal T = lData[EulerTerm::T];
  
  // this is inconsistent with Prabhu's linearization
  //  const CFreal rho = lData[EulerTerm::RHO];
  //const CFreal p = lData[EulerTerm::P];
  //const CFreal cvTr = eData->dEdT;
  //const CFreal beta = p/(rho*T)/cvTr;
  
  const CFuint firstSpecies = getModel()->getFirstScalarVar(0);
  CFreal sumAlphaYs = 0.0;
  
  cf_assert(_ys.size() == nbSpecies);
  cf_assert(_alpha.size() == nbSpecies);
  //  cf_assert(p > 0.0);
  cf_assert(T > 0.0);
  //cf_assert(rho > 0.0);

  CFreal numBeta = 0.;
  CFreal denBeta = 0.;
  for (CFuint i = 0; i < nbSpecies; ++i) {
    const CFreal sigmai = lData[firstSpecies + i]/_mmasses[i];
    numBeta += sigmai;
    denBeta += sigmai*_fcoeff[i];
  }
  
  const CFreal beta = numBeta/denBeta;
  for (CFuint is = 0; is < nbSpecies; ++is) {
    _ys[is] = lData[firstSpecies + is];
    
    if (_ys[is] > 1.1) {
      cout << "_ys > 1.1 = " << _ys << endl;
      // cf_assert(_ys[is] <= 1.1);
    }
    _alpha[is] = _RiGas[is]*T - beta*((eData->energyTr)[is]);
    sumAlphaYs += _alpha[is]*_ys[is];
  }

  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal u = lData[EulerTerm::VX];
  const CFreal v = lData[EulerTerm::VY];

  const CFreal rho = lData[EulerTerm::RHO];//My addition

  const CFreal V2 = lData[EulerTerm::V]*lData[EulerTerm::V];
  const CFreal q = 0.5*V2;
  const CFreal V = v*nx -u*ny;
  const CFreal H = lData[EulerTerm::H];
  const CFuint firstTv = getModel()->getFirstScalarVar(1);
  const CFreal ev = lData[firstTv];
  //  const CFreal a2 = (1. + beta)*p/rho;
  // const CFreal a = sqrt(a2);
  const CFreal a = lData[EulerTerm::A];
  const CFreal a2 = lData[EulerTerm::A]*lData[EulerTerm::A];
  const CFreal U = u*nx + v*ny;
  
  const CFuint nbSpPlus1 = nbSpecies+1;
  const CFuint nbSpPlus2 = nbSpecies+2;
  const CFuint nbSpPlus3 = nbSpecies+3;

  // rightEv = 0.0;
  // leftEv = 0.0;
  for (CFuint is = 0; is < nbSpecies; ++is) {
    rightEv(is, is) = -rho;
    rightEv(is,nbSpPlus1) = 0.0;//Is the matrix rightEv initialized as zero. Otherwise I would need to define rightEv(is,nbSp)=0.0
    rightEv(is,nbSpPlus2) = 0.0;

    rightEv(nbSpecies,is) = 0.0;
    rightEv(nbSpPlus1,is) = 0.0;
    rightEv(nbSpPlus2,is) = 0.0;
    rightEv(nbSpPlus3,is) = ev;//This is my addition
  }

  rightEv(nbSpecies,nbSpecies) = -ny;
  rightEv(nbSpecies,nbSpPlus1) = 0.5*nx/a;
  rightEv(nbSpecies,nbSpPlus2) = -0.5*nx/a;
  //Is the matrix rightEv initialized as zero. Otherwise I would need to define rightEv(nbSpecies,nbSpPlus3)=0.0
  rightEv(nbSpPlus1,nbSpecies) = nx;
  rightEv(nbSpPlus1,nbSpPlus1) = 0.5*ny/a;; // 0.5(v+a*ny) ???
  rightEv(nbSpPlus1,nbSpPlus2) = -0.5*ny/a;; // 0.5(v+a*ny) ???

  rightEv(nbSpPlus2,nbSpecies) = 0.0;
  rightEv(nbSpPlus2,nbSpPlus1) = 0.5/a;
  rightEv(nbSpPlus2,nbSpPlus2) = 0.5/a;
  rightEv(nbSpPlus2,nbSpPlus3) = 0.0;

  rightEv(nbSpPlus3,nbSpPlus1) = 0.5*ev;
  rightEv(nbSpPlus3,nbSpPlus2) = 0.5*ev;
  
  rightEv(nbSpPlus3,nbSpPlus3) = -1.0;

  rightEv /= rho;//rightEv checked on 2009/3/30.

  //  cout << "R = "<< endl;
  //   cout << rightEv << endl <<endl;

  const CFreal bq = beta*q;
  for (CFuint is = 0; is < nbSpecies; ++is) {
    for (CFuint js = 0; js < nbSpecies; ++js) {
      const CFreal minusdelta = (js != is) ? 0.0 : -1.0;
      leftEv(is,js) = minusdelta;
    }//This could be simplified.

    leftEv(is,nbSpecies) = 0.0;
    leftEv(is,nbSpPlus1) = 0.0;
    leftEv(is,nbSpPlus2) = 0.0;
    leftEv(is,nbSpPlus3) = 0.0;

    const CFreal alphaPlusBetaq = _alpha[is] + bq;
    leftEv(nbSpecies,is) = 0.0;
    leftEv(nbSpPlus1,is) = 0.0;
    leftEv(nbSpPlus2,is) = 0.0;
    leftEv(nbSpPlus3,is) = -ev;
  }

  leftEv(nbSpecies,nbSpecies) = -rho*ny;
  leftEv(nbSpecies,nbSpPlus1) = rho*nx;

  leftEv(nbSpPlus1,nbSpecies) = rho*a*nx; //a*nx - beta*u;
  leftEv(nbSpPlus1,nbSpPlus1) = rho*a*ny; //a*ny - beta*v;
  leftEv(nbSpPlus1,nbSpPlus2) = rho*a;
  leftEv(nbSpPlus1,nbSpPlus3) = 0.0;

  leftEv(nbSpPlus2,nbSpecies) = -rho*a*nx;
  leftEv(nbSpPlus2,nbSpPlus1) = -rho*a*ny;
  leftEv(nbSpPlus2,nbSpPlus2) = rho*a;
  leftEv(nbSpPlus2,nbSpPlus3) = 0.0;

  leftEv(nbSpPlus3,nbSpecies) = 0.0;
  leftEv(nbSpPlus3,nbSpPlus1) = 0.0;
  leftEv(nbSpPlus3,nbSpPlus2) = rho*ev*a;
  leftEv(nbSpPlus3,nbSpPlus3) = -rho;//leftEv checked on 2009/3/30.

  eValues = U;
  eValues[nbSpPlus1] += a;
  eValues[nbSpPlus2] -= a;//eValues checked on 2009/3/30.
}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler2DNEQSymm::getBlockSeparator() const
{
  return Framework::PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm::splitJacobian(RealMatrix& jacobPlus,
                                RealMatrix& jacobMin,
                                RealVector& eValues,
                                const RealVector& normal)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint nbTv = getModel()->getNbScalarVars(1);
  cf_assert(nbTv == 1);
  
  const RealVector& lData = getModel()->getPhysicalData();
  SafePtr<PhysicalChemicalLibrary::ExtraData> eData = _library->getExtraData();
    
  const CFreal T = lData[EulerTerm::T];
  const CFreal rho = lData[EulerTerm::RHO];
  const CFreal p = lData[EulerTerm::P];
  const CFreal cvTr = eData->dEdT;
  const CFreal beta = p/(rho*T)/cvTr;
  const CFuint firstSpecies = getModel()->getFirstScalarVar(0);
  CFreal sumAlphaYs = 0.0;
  
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
    
  for (CFuint is = 0; is < nbSpecies; ++is) {
    _ys[is] = lData[firstSpecies + is];
    
    /*if (_ys[is] > 1.0) {
      cout << "_ys = " << _ys << endl;
      cout << "rho = " << rho << endl;
      cout << "T   = " << T << endl;
      cout << "p   = " << p << endl;
      cout << "cvTr= " << cvTr << endl << endl;
      //  cf_assert(_ys[is] <= 1.0);
    }*/
    
    _alpha[is] = _RiGas[is]*T - beta*(eData->energyTr)[is];
    sumAlphaYs += _alpha[is]*_ys[is];
  }  
  
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal u = lData[EulerTerm::VX];
  const CFreal v = lData[EulerTerm::VY];
  const CFreal V2 = lData[EulerTerm::V]*lData[EulerTerm::V];
  const CFreal q = 0.5*V2;
  const CFreal V = v*nx -u*ny;
  const CFreal H = lData[EulerTerm::H];
  const CFuint firstTv = getModel()->getFirstScalarVar(1);
  const CFreal ev = lData[firstTv];
  const CFreal a2 = (1. + beta)*p/rho;
  const CFreal a = sqrt(a2);
  const CFreal U = u*nx + v*ny;
  const CFuint nbSpPlus1 = nbSpecies+1; 
  const CFuint nbSpPlus2 = nbSpecies+2; 
  const CFuint nbSpPlus3 = nbSpecies+3;
  
  // _rightEv = 0.0;
  // _leftEv = 0.0;
  for (CFuint is = 0; is < nbSpecies; ++is) {
    _rightEv(is, is) = -rho;
    _rightEv(is,nbSpPlus1) = 0.0;
    _rightEv(is,nbSpPlus2) = 0.0;
    //Is the matrix rightEv initialized to zero? Otherwise I would need to define rightEv(is,nbSp)=0.0
    _rightEv(nbSpecies,is) = 0.0;
    _rightEv(nbSpPlus1,is) = 0.0;
    _rightEv(nbSpPlus2,is) = 0.0;
    _rightEv(nbSpPlus3,is) = ev;//My addition.
  }  
  
  _rightEv(nbSpecies,nbSpecies) = -ny;
  _rightEv(nbSpecies,nbSpPlus1) = 0.5*nx/a;
  _rightEv(nbSpecies,nbSpPlus2) = -0.5*nx/a;
  
  _rightEv(nbSpPlus1,nbSpecies) = nx;
  _rightEv(nbSpPlus1,nbSpPlus1) = 0.5*ny/a; // 0.5(v+a*ny) ???
  _rightEv(nbSpPlus1,nbSpPlus2) = -0.5*ny/a; // 0.5(v+a*ny) ???
  
  _rightEv(nbSpPlus2,nbSpecies) = 0.0;
  _rightEv(nbSpPlus2,nbSpPlus1) = 0.5/a;
  _rightEv(nbSpPlus2,nbSpPlus2) = 0.5/a;
  _rightEv(nbSpPlus2,nbSpPlus3) = 0.0;
  
  _rightEv(nbSpPlus3,nbSpPlus1) = 0.5*ev;
  _rightEv(nbSpPlus3,nbSpPlus2) = 0.5*ev;
  _rightEv(nbSpPlus3,nbSpPlus3) = -1.0;
  
  _rightEv /= rho;//_rightEv checked on 2009/3/30.
  
 //  cout << "R = "<< endl;
//   cout << _rightEv << endl <<endl;
  
  const CFreal bq = beta*q;
  for (CFuint is = 0; is < nbSpecies; ++is) {
    for (CFuint js = 0; js < nbSpecies; ++js) {
      const CFreal minusdelta = (js != is) ? 0.0 : -1.0; 
      _leftEv(is,js) = minusdelta;
    }//This could be simplified.
    
    _leftEv(is,nbSpecies) = 0.0;
    _leftEv(is,nbSpPlus1) = 0.0;
    _leftEv(is,nbSpPlus2) = 0.0;
    _leftEv(is,nbSpPlus3) = 0.0;
    
    //const CFreal alphaPlusBetaq = _alpha[is] + bq;
    _leftEv(nbSpecies,is) = 0.0;
    _leftEv(nbSpPlus1,is) = 0.0;
    _leftEv(nbSpPlus2,is) = 0.0;
    _leftEv(nbSpPlus3,is) = -ev;
  }  
  
  _leftEv(nbSpecies,nbSpecies) = -rho*ny;
  _leftEv(nbSpecies,nbSpPlus1) = rho*nx;
  
  _leftEv(nbSpPlus1,nbSpecies) = rho*a*nx; //a*nx - beta*u;
  _leftEv(nbSpPlus1,nbSpPlus1) = rho*a*ny; //a*ny - beta*v;
  _leftEv(nbSpPlus1,nbSpPlus2) = rho*a;
  _leftEv(nbSpPlus1,nbSpPlus3) = 0.0;
  
  _leftEv(nbSpPlus2,nbSpecies) = -rho*a*nx;
  _leftEv(nbSpPlus2,nbSpPlus1) = -rho*a*ny;
  _leftEv(nbSpPlus2,nbSpPlus2) = rho*a;
  _leftEv(nbSpPlus2,nbSpPlus3) = 0.0;
  
  _leftEv(nbSpPlus3,nbSpecies) = 0.0; 
  _leftEv(nbSpPlus3,nbSpPlus1) = 0.0;
  _leftEv(nbSpPlus3,nbSpPlus2) = rho*ev*a;
  _leftEv(nbSpPlus3,nbSpPlus3) = -rho;//_leftEv checked on 2009/3/30.
  
  //Comment back from here to
  //cout << "L = "<< endl;
  // cout << _leftEv << endl <<endl;
  
 // RealMatrix mat(_rightEv.nbRows(), _leftEv.nbRows());
  //mat = _rightEv*_leftEv;
  //cout <<"R*L" << endl;
  //cout << mat << endl << endl;
  
  // mat = _leftEv*_rightEv;
  // cout <<"L*R" << endl;
  // cout << mat << endl << endl;
  //here
  
  eValues = U;
  eValues[nbSpPlus1] += a;
  eValues[nbSpPlus2] -= a;
  
  //Comment back from here to:
  /*
  RealMatrix mat(_rightEv.nbRows(), _leftEv.nbRows());
  //const CFreal bq = beta*q;
  
  for (CFuint is = 0; is < nbSpecies; ++is) {
    for (CFuint js = 0; js < nbSpecies; ++js) {
      const CFreal Udelta = (js != is) ? 0.0 : U; 
      mat(is,js) = Udelta;
    }
    mat(is,nbSpecies) = 0.0;
    mat(is,nbSpPlus1) = 0.0;
    mat(is,nbSpPlus2) = 0.0;
    mat(is,nbSpPlus3) = 0.0;
    
    const CFreal alphaPlusBeta = _alpha[is] + bq;
    mat(nbSpecies,is) = 0.0;
    mat(nbSpPlus1,is) = 0.0;
    mat(nbSpPlus2,is) = 0.0;
    mat(nbSpPlus3,is) = 0.0;
  }
  
  mat(nbSpecies,nbSpecies) = U;
  mat(nbSpecies,nbSpPlus1) = 0.0;
  mat(nbSpecies,nbSpPlus2) = a*nx;
  mat(nbSpecies,nbSpPlus3) = 0.0;
  
  mat(nbSpPlus1,nbSpecies) = 0.0;
  mat(nbSpPlus1,nbSpPlus1) = U;
  mat(nbSpPlus1,nbSpPlus2) = a*ny;
  mat(nbSpPlus1,nbSpPlus3) = 0.0;
  
  mat(nbSpPlus2,nbSpecies) = a*nx;
  mat(nbSpPlus2,nbSpPlus1) = a*ny;
  mat(nbSpPlus2,nbSpPlus2) = U;
  mat(nbSpPlus2,nbSpPlus3) = 0.0;
  
  mat(nbSpPlus3,nbSpecies) = ev*a2*nx;
  mat(nbSpPlus3,nbSpPlus1) = ev*a2*ny;
  mat(nbSpPlus3,nbSpPlus2) = 0.0;
  mat(nbSpPlus3,nbSpPlus3) = U;//mat checked on 2009/3/30.

   RealMatrix m1(9,9);
   RealMatrix m2(9,9);
  
   m1 = mat*_rightEv;
   m2 = _leftEv*m1;
  
   cout << "lambda M = " << endl;
   cout << m2 << endl;
   cout << "lambda = " << eValues << endl << endl; 
   //to here: the matrix m2 is diagonal with entries the eigenvalues of the convective jacobian; checked on 2009/3/30.
   */

// jacobian matrix
//   mat = 0.0;
//   for (CFuint is = 0; is < nbSpecies; ++is) {
//     for (CFuint js = 0; js < nbSpecies; ++js) {
//       const CFreal delta = (js != is) ? 0.0 : 1.0; 
//       mat(is,js) = (delta - _ys[is])*U;
//     }
//     mat(is,nbSpecies) = _ys[is]*nx;
//     mat(is,nbSpPlus1) = _ys[is]*ny;
//     mat(is,nbSpPlus2) = 0.0;
//     mat(is,nbSpPlus3) = 0.0;
    
//     const CFreal alphaPlusBeta = _alpha[is] + bq;
//     mat(nbSpecies,is) = alphaPlusBeta*nx -U*u;
//     mat(nbSpPlus1,is) = alphaPlusBeta*ny -U*v;
//     mat(nbSpPlus2,is) = (alphaPlusBeta - H)*U;
//     mat(nbSpPlus3,is) = -U*ev;
//   }
  
//   mat(nbSpecies,nbSpecies) = (1.-beta)*u*nx+U;
//   mat(nbSpecies,nbSpPlus1) = (1.-beta)*v*nx-V;
//   mat(nbSpecies,nbSpPlus2) = beta*nx;
//   mat(nbSpecies,nbSpPlus3) = -beta*nx;
  
//   mat(nbSpPlus1,nbSpecies) = (1.-beta)*u*ny+V;
//   mat(nbSpPlus1,nbSpPlus1) = (1.-beta)*v*ny+U;
//   mat(nbSpPlus1,nbSpPlus2) = beta*ny;
//   mat(nbSpPlus1,nbSpPlus3) = -beta*ny;
  
//   mat(nbSpPlus2,nbSpecies) = H*nx - beta*U*u;
//   mat(nbSpPlus2,nbSpPlus1) = H*ny - beta*U*v;
//   mat(nbSpPlus2,nbSpPlus2) = (1.+beta)*U;
//   mat(nbSpPlus2,nbSpPlus3) = -beta*U;
  
//   mat(nbSpPlus3,nbSpecies) = ev*nx;
//   mat(nbSpPlus3,nbSpPlus1) = ev*ny;
//   mat(nbSpPlus3,nbSpPlus2) = 0.0;
//   mat(nbSpPlus3,nbSpPlus3) = U;
  
//   RealMatrix m1(9,9);
//   RealMatrix m2(9,9);
  
//   m1 = mat*_rightEv;
//   m2 = _leftEv*m1;
  
//   cout << "lambda M = " << endl;
//   cout << m2 << endl;
//   cout << "lambda = " << eValues << endl << endl;
  //here


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
  CFLogDebugMax( "RightEigenvectors @Euler2DNEQSymm::splitJacobian" << "\n"
		 << _rightEv << "\n");
  CFLogDebugMax( "LeftEigenvectors @Euler2DNEQSymm::splitJacobian" << "\n"
		 << _leftEv << "\n");
  CFLogDebugMax( "EigenValues @Euler2DNEQSymm::splitJacobian" << "\n"
		 << eValues << "\n" << "\n"); 
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm::computePhysicalData(State& state,
				     RealVector& data)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  CFreal rho = 0.0;
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    rho += state[ie];
  }
  
  data[EulerTerm::RHO] = rho;
  const CFreal ovRho = 1./rho;
  
  // set the species mass fractions
  const CFuint firstSpecies = getModel()->getFirstScalarVar(0);
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    _ys[ie] = state[ie]*ovRho;
    data[firstSpecies + ie] = _ys[ie];
  }

  // U  and V Velocity Average
  data[EulerTerm::VX] = state[nbSpecies]*ovRho;
  data[EulerTerm::VY] = state[nbSpecies+1]*ovRho;
  
  const CFreal V2 = data[EulerTerm::VX]*data[EulerTerm::VX] +
    data[EulerTerm::VY]*data[EulerTerm::VY];
  data[EulerTerm::V] = sqrt(V2);
  
  const CFuint startTv = getModel()->getFirstScalarVar(1);
  data[startTv]= state[nbSpecies+3]*ovRho; //ev

  data[EulerTerm::E]  = state[nbSpecies+2]*ovRho;
  const CFreal Rgas = _library->getRgas();
  SafePtr<PhysicalChemicalLibrary::ExtraData> eData = _library->getExtraData();
  
  CFreal denom = 0.;
  CFreal form  = 0.;
  CFreal riovermi  = 0.;
  for (CFuint i = 0; i < nbSpecies; ++i) {
    riovermi += state[i]/_mmasses[i];
    const CFreal yOvM = _ys[i]/_mmasses[i];
    denom += yOvM*((Rgas*_fcoeff[i]));
    form += _ys[i]*eData->enthalpyForm[i];
  }
  
  data[EulerTerm::T] = (data[EulerTerm::E] - data[startTv] - form-0.5*V2)/denom;

  const CFreal P = data[EulerTerm::T]*Rgas*riovermi;
  data[EulerTerm::P] = P;
  data[EulerTerm::H] = data[EulerTerm::E] + P*ovRho;
 
  //Speed of Sound  Average
  CFreal numBeta = 0.;
  CFreal denBeta = 0.;
  for (CFuint i = 0; i < nbSpecies; ++i) {
    const CFreal sigmai = _ys[i]/_mmasses[i];
    numBeta += sigmai;
    denBeta += sigmai*_fcoeff[i];
  }

  const CFreal beta = numBeta/denBeta;
  const CFreal RT = Rgas*data[EulerTerm::T];

  CFreal aiyi = 0.0;
  for (CFuint i = 0; i < nbSpecies; ++i) {
    aiyi += (_ys[i]/_mmasses[i])*(RT - beta*(_fcoeff[i]*RT + _mmasses[i]*eData->enthalpyForm[i]));
  }
  
  data[EulerTerm::A] = std::sqrt(aiyi + beta*(data[EulerTerm::H] - 0.5*V2 - data[startTv]));
  // cout << "data = " << data<< endl;
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm::computeStateFromPhysicalData(const RealVector& data,
					 State& state)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DNEQSymm::computeStateFromPhysicalData()");
  
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

void Euler2DNEQSymm::computePerturbedStatesData
(const vector<State*>& states,
 const CFuint nbStatesInVec,
 const CFuint iVar)
{
  throw Common::NotImplementedException (FromHere(),"Euler2DNEQSymm::computePerturbedStatesData()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQSymm::computeProjectedJacobian(const RealVector& normal, RealMatrix& jacob)
  {
  //Jacobian in Qvars.
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  const CFuint nbTv = getModel()->getNbScalarVars(1);
  cf_assert(nbTv == 1);
  
  const RealVector& lData = getModel()->getPhysicalData();
  SafePtr<PhysicalChemicalLibrary::ExtraData> eData = _library->getExtraData();
  
  const CFreal T = lData[EulerTerm::T];
  const CFreal rho = lData[EulerTerm::RHO];
  const CFreal p = lData[EulerTerm::P];
  const CFreal cvTr = eData->dEdT;
  const CFreal beta = p/(rho*T)/cvTr;
  const CFuint firstSpecies = getModel()->getFirstScalarVar(0);
  CFreal sumAlphaYs = 0.0;
  
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
  
  for (CFuint is = 0; is < nbSpecies; ++is) {
    _ys[is] = lData[firstSpecies + is];
    cf_assert(_ys[is] < 1.1);
    
    _alpha[is] = _RiGas[is]*T - beta*(eData->energyTr)[is];
    sumAlphaYs += _alpha[is]*_ys[is];
  }  
  
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal u = lData[EulerTerm::VX];
  const CFreal v = lData[EulerTerm::VY];
  const CFreal V2 = lData[EulerTerm::V]*lData[EulerTerm::V];
  const CFreal q = 0.5*V2;
  const CFreal V = v*nx -u*ny;
  const CFreal H = lData[EulerTerm::H];
  const CFuint firstTv = getModel()->getFirstScalarVar(1);
  const CFreal ev = lData[firstTv];
  const CFreal a2 = (1. + beta)*p/rho;

  const CFreal a = sqrt(a2);//JGM's addition.

  const CFreal U = u*nx + v*ny;
  const CFuint nbSpPlus1 = nbSpecies+1; 
  const CFuint nbSpPlus2 = nbSpecies+2; 
  const CFuint nbSpPlus3 = nbSpecies+3;
  const CFreal bq = beta*q;
  
  for (CFuint is = 0; is < nbSpecies; ++is) {
    for (CFuint js = 0; js < nbSpecies; ++js) {
      const CFreal Udelta = (js != is) ? 0.0 : U; 
      jacob(is,js) = Udelta;
    }
    jacob(is,nbSpecies) = 0.0;
    jacob(is,nbSpPlus1) = 0.0;
    jacob(is,nbSpPlus2) = 0.0;
    jacob(is,nbSpPlus3) = 0.0;
    
    const CFreal alphaPlusBeta = _alpha[is] + bq;
    jacob(nbSpecies,is) = 0.0;
    jacob(nbSpPlus1,is) = 0.0;
    jacob(nbSpPlus2,is) = 0.0;
    jacob(nbSpPlus3,is) = 0.0;
  }
  
  jacob(nbSpecies,nbSpecies) = U;
  jacob(nbSpecies,nbSpPlus1) = 0.0;
  jacob(nbSpecies,nbSpPlus2) = a*nx;
  jacob(nbSpecies,nbSpPlus3) = 0.0;
  
  jacob(nbSpPlus1,nbSpecies) = 0.0;
  jacob(nbSpPlus1,nbSpPlus1) = U;
  jacob(nbSpPlus1,nbSpPlus2) = a*ny;
  jacob(nbSpPlus1,nbSpPlus3) = 0.0;
  
  jacob(nbSpPlus2,nbSpecies) = a*nx;
  jacob(nbSpPlus2,nbSpPlus1) = a*ny;
  jacob(nbSpPlus2,nbSpPlus2) = U;
  jacob(nbSpPlus2,nbSpPlus3) = 0.0;
  
  jacob(nbSpPlus3,nbSpecies) = ev*a2*nx;
  jacob(nbSpPlus3,nbSpPlus1) = ev*a2*ny;
  jacob(nbSpPlus3,nbSpPlus2) = 0.0;
  jacob(nbSpPlus3,nbSpPlus3) = U;//jacob checked on 2009/3/30.
//Computing jacobian ended.
  
  //Now we need dU/dQ and dQ/dU.
  // First define dU/dQ, that is, SymmToCons:
  //cf_assert(_model.isNotNull());
  //const RealVector& lData = _model->getPhysicalData();
  
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  //Common::SafePtr<PhysicalChemicalLibrary::ExtraData> eData = library->getExtraData();
  
  //const CFuint nbSpecies = _model->getNbScalarVars(0);
  //const CFreal rho = linearData[EulerTerm::RHO];
  const CFreal ovrho = 1./rho;
  //const CFreal u = linearData[EulerTerm::VX];
  //const CFreal v = linearData[EulerTerm::VY];
  //JGM's additions:
  //const CFreal a = linearData[EulerTerm::A];
  //const CFreal a2 = linearData[EulerTerm::A]*linearData[EulerTerm::A];
  const CFreal ova = 1./a;
  const CFreal ova2 = 1./a2;
  ////////////////////////
  const CFuint uID  = nbSpecies;
  const CFuint vID  = nbSpecies+1;
  const CFuint eID  = nbSpecies+2;
  const CFuint evID = nbSpecies+3;
  const CFreal eT = eData->dEdT; // check this !!!!
  // DevDTv is stored in PhysicalChemicalLibrary during linearization
  const CFreal evTv = eData->dEvTv;
  //Would be redefinition: const CFuint firstTv = _model->getFirstScalarVar(1);
  cf_assert(_model->getNbScalarVars(1) == 1);
  

  //JGM's additions:
    
  //Would be redefinition: const CFuint firstSpecies = _model->getFirstScalarVar(0);
  //Would be redefinition: CFreal sumAlphaYs = 0.0;

  /*_ys.resize(nbSpecies);
  _alpha.resize(nbSpecies);
  _RiGas.resize(nbSpecies);
  _mmasses.resize(nbSpecies);
  _fcoeff.resize(nbSpecies);*/
  
  //The previous vectors are initialized:

  cf_assert(_ys.size() == nbSpecies);
  cf_assert(_alpha.size() == nbSpecies);
  //  cf_assert(p > 0.0);
  //cf_assert(T > 0.0);
  //cf_assert(rho > 0.0);
  
  library->setRiGas(_RiGas);//Think on commenting this.

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

  /*
  CFreal numBeta = 0.;
  CFreal denBeta = 0.;
  for (CFuint i = 0; i < nbSpecies; ++i) {
    const CFreal sigmai = linearData[firstSpecies + i]/_mmasses[i];
    numBeta += sigmai;
    denBeta += sigmai*_fcoeff[i];
  }*/

  //const CFreal T = linearData[EulerTerm::T];

  //const CFreal beta = numBeta/denBeta;

  /*for (CFuint is = 0; is < nbSpecies; ++is) {
    _ys[is] = linearData[firstSpecies + is];

    if (_ys[is] > 1.1) {
      cout << "_ys > 1.1 = " << _ys << endl;
      // cf_assert(_ys[is] <= 1.1);
    }
  
    _alpha[is] = _RiGas[is]*T - beta*((eData->energyTr)[is]);
    sumAlphaYs += _alpha[is]*_ys[is];
  }*/
  // Vectors initialization ends here.
  
  /*const CFreal V2 = linearData[EulerTerm::V]*linearData[EulerTerm::V];
  const CFreal ev = linearData[firstTv];*/

  /////////////////////////////////////////////////////////
  RealMatrix dUdQ(9,9);
  for (CFuint countRows=0;countRows<9; ++countRows){
    for (CFuint countColumns=0;countColumns<9; ++countColumns){
      dUdQ(countRows,countColumns)=0.0;
    }
  }

  for (CFuint is = 0; is < nbSpecies; ++is) {
    dUdQ(is,is) = -ova2;
    dUdQ(uID,is) = -ova2*u;
    dUdQ(vID,is) = -ova2*v;
    dUdQ(eID,is) = -ova2*(ev+V2)+ova2*_alpha[is]/beta;
    dUdQ(evID,is) = -ova2*ev;
  }
  
for (CFuint is = 0; is < nbSpecies; ++is) {
  dUdQ(is,eID) = ova*rho*_ys[is];
 }

  dUdQ(uID,uID) = rho;
  dUdQ(vID,vID) = rho;
  
  dUdQ(uID,eID) = rho*ova*u;
  dUdQ(vID,eID) = rho*ova*v;

  dUdQ(eID,uID)  = rho*u;
  dUdQ(eID,vID)  = rho*v;
  dUdQ(eID,eID)  = ova*rho*(V2+2*ev)+rho*a/beta-ova*rho*sumAlphaYs/beta;
  dUdQ(eID,evID) = -rho*ova2;
  
  dUdQ(uID,evID) = 0.0;//Should these 
  dUdQ(vID,evID) = 0.0;//two terms be not defined????
  dUdQ(eID,evID) = 2*rho*ev*ova;
  dUdQ(evID,evID) = -rho*ova2;//_transMatrix checked on 2009/4/2. Error in position (eID,eID) fixed.

   //Remove this:
   /*cout << "Enter into SymmToCons: ";
   cout << "dU/dQ = " << endl;
   cout << dUdQ << endl;*/
  //End of dU/dQ

  //Then define dQ/dU, conversely ConsToSymm:

  RealMatrix dQdU(9,9);
  for (CFuint countRows=0;countRows<9; ++countRows){
    for (CFuint countColumns=0;countColumns<9; ++countColumns){
      dQdU(countRows,countColumns)=0.0;
    }
  }

  for (CFuint is = 0; is < nbSpecies; ++is) {
    for (CFuint js = 0; js < nbSpecies; ++js) {
      dQdU(is,js) = _ys[is]*_alpha[js];
      if (is==js)
	dQdU(is,js) =dQdU(is,is)-a2;
    }

    dQdU(uID,is) = -ovrho*u;
    dQdU(vID,is) = -ovrho*v;
    dQdU(eID,is) = ovrho*ova*_alpha[is];
    dQdU(evID,is) = ovrho*ev*(_alpha[is]+a2);
  }
  
    for (CFuint is = 0; is < nbSpecies; ++is) {
    dQdU(is,uID) = -beta*_ys[is]*u;
    dQdU(is,vID) = -beta*_ys[is]*v;
    dQdU(is,eID) = beta*_ys[is];
    dQdU(is,evID) = -beta*_ys[is];//We are assuming phi=-beta, that is, no ionization of the plasma.
  }

  dQdU(uID,uID) = ovrho;
  dQdU(vID,vID) = ovrho;
  
  dQdU(eID,uID) = -ovrho*ova*beta*u;
  dQdU(eID,vID) = -ovrho*ova*beta*v;
  dQdU(eID,eID)  = ovrho*ova*beta;
  dQdU(eID,evID) = -ovrho*ova*beta;
  
  dQdU(uID,evID) = -ovrho*beta*ev*u;
  dQdU(vID,evID) = -ovrho*beta*ev*v;
  dQdU(eID,evID) = ovrho*beta*ev;      
  dQdU(evID,evID) = ovrho*(-beta*ev-a2);//dQdU checked on 2009/3/30.
  //End of dQ/dU
   
   /*
   //Remove this:
   cout << "Enter into ConsToSymm dQ/dU = " << endl;
   cout << dQdU << endl;
   */
//Define jacobU:
 //jacobian matrix in terms of U.
   RealMatrix jacobU;
   for (CFuint countRows=0;countRows<9; ++countRows){
    for (CFuint countColumns=0;countColumns<9; ++countColumns){
      jacobU(countRows,countColumns)=0.0;
    }
  }

   for (CFuint is = 0; is < nbSpecies; ++is) {
     for (CFuint js = 0; js < nbSpecies; ++js) {
       const CFreal delta = (js != is) ? 0.0 : 1.0; 
       jacobU(is,js) = (delta - _ys[is])*U;
     }
     jacobU(is,nbSpecies) = _ys[is]*nx;
     jacobU(is,nbSpPlus1) = _ys[is]*ny;
     jacobU(is,nbSpPlus2) = 0.0;
     jacobU(is,nbSpPlus3) = 0.0;
    
     //const CFreal alphaPlusBeta = _alpha[is] + bq;
     jacobU(nbSpecies,is) = _alpha[is]*nx -U*u;
     jacobU(nbSpPlus1,is) = _alpha[is]*ny -U*v;
     jacobU(nbSpPlus2,is) = (_alpha[is] - H)*U;
     jacobU(nbSpPlus3,is) = -U*ev;
   }
  
   jacobU(nbSpecies,nbSpecies) = (1.-beta)*u*nx+U;
   jacobU(nbSpecies,nbSpPlus1) = v*nx-beta*ny*u;
   jacobU(nbSpecies,nbSpPlus2) = beta*nx;//Check the sign.
   jacobU(nbSpecies,nbSpPlus3) = -beta*nx;
  
   jacobU(nbSpPlus1,nbSpecies) = v*nx-beta*ny*u;
   jacobU(nbSpPlus1,nbSpPlus1) = (1.-beta)*v*ny+U;
   jacobU(nbSpPlus1,nbSpPlus2) = beta*ny;//Check the sign.
   jacobU(nbSpPlus1,nbSpPlus3) = -beta*ny;
  
   jacobU(nbSpPlus2,nbSpecies) = H*nx - beta*U*u;
   jacobU(nbSpPlus2,nbSpPlus1) = H*ny - beta*U*v;
   jacobU(nbSpPlus2,nbSpPlus2) = (1.+beta)*U;
   jacobU(nbSpPlus2,nbSpPlus3) = -beta*U;
  
   jacobU(nbSpPlus3,nbSpecies) = ev*nx;
   jacobU(nbSpPlus3,nbSpPlus1) = ev*ny;
   jacobU(nbSpPlus3,nbSpPlus2) = 0.0;
   jacobU(nbSpPlus3,nbSpPlus3) = U;
//End of jacobU

   RealMatrix m1(9,9);
   RealMatrix m2(9,9);
  
   m1 = jacob*dQdU;
   m2 = dUdQ*m1;
   cout << "dUdQ*jacob*dQdU = "<< endl;
   cout << m2 << endl <<endl;
  
   cout << "jacobU = "<< endl;
   cout << jacobU << endl <<endl;

}


//////////////////////////////////////////////////////////////////////////////

    } // end of namespace NEQ

  } // end of namespace Physics

} // end of namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
