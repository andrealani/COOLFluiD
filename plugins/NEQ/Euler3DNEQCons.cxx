#include "NEQ.hh"
#include "Euler3DNEQCons.hh"
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

Environment::ObjectProvider<Euler3DNEQCons, ConvectiveVarSet, NEQModule, 1>
eulerNEQ3DConsProvider("Euler3DNEQCons");

//////////////////////////////////////////////////////////////////////////////

Euler3DNEQCons::Euler3DNEQCons(Common::SafePtr<BaseTerm> term) :
  MultiScalarVarSet<Euler3DVarSet>(term),
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
  
  vector<std::string> names(nbSpecies + 4 + nbTv);
  for (CFuint ie = 0; ie < nbSpecies; ++ie) {
    names[ie] = "rho" + StringOps::to_str(ie);
  }
  names[nbSpecies]     = "rhoU";
  names[nbSpecies + 1] = "rhoV";
  names[nbSpecies + 2] = "rhoW";
  names[nbSpecies + 3] = "rhoE";
  
  const CFuint startTv = nbSpecies + 4; 
  for (CFuint ie = 0; ie < nbTv; ++ie) {
    names[startTv + ie] = "rhoEv" + StringOps::to_str(ie);
  }

  setVarNames(names);
}

//////////////////////////////////////////////////////////////////////////////

Euler3DNEQCons::~Euler3DNEQCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQCons::computeEigenValuesVectors(RealMatrix& rightEv,
                                           RealMatrix& leftEv,
                                           RealVector& eValues,
                                           const RealVector& normal)
{
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  cf_assert(getModel()->getNbScalarVars(1) == 1);

  const RealVector& lData = getModel()->getPhysicalData();
  SafePtr<PhysicalChemicalLibrary::ExtraData> eData = _library->getExtraData();

  const CFreal T = lData[EulerTerm::T];

  // this is inconsistent with Prabhu's linearization
  //  const CFreal rho = lData[EulerTerm::RHO];
  //const CFreal p = lData[EulerTerm::P];
  //const CFreal cvTr = eData->dEdT;
  //const CFreal beta = p/(rho*T)/cvTr;

  const CFuint firstSpecies = getModel()->getFirstScalarVar(0);
  
  cf_assert(_ys.size() == nbSpecies);
  cf_assert(_alpha.size() == nbSpecies);
  //  cf_assert(p > 0.0);
  cf_assert(T > 0.0);
  //cf_assert(rho > 0.0);

  CFreal numBeta = 0.;
  CFreal denBeta = 0.;
  const CFuint start = (_library->presenceElectron()) ? 1 : 0;
  for (CFuint i = start; i < nbSpecies; ++i) {
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
    if (!_library->presenceElectron()) {
      Tq = T;
    }
    else {
      if (is == 0) {
	Tq = 0.0; // get vibrational temperature
      }
      else{
	Tq = T;
      }
    }
    
    _alpha[is] = _RiGas[is]*Tq - beta*(eData->energyTr)[is]; // change here to include phi
  }
  
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal nz = normal[ZZ];
  const CFreal u = lData[EulerTerm::VX];
  const CFreal v = lData[EulerTerm::VY];
  const CFreal w = lData[EulerTerm::VZ];

  const CFreal V2 = lData[EulerTerm::V]*lData[EulerTerm::V];
  const CFreal q = 0.5*V2;

  // defining wall tangentional components

  // arbitrary unit vector (not parallel to unit normal)
  const CFreal vx = nz;
  const CFreal vy = nx;
  const CFreal vz = ny;

  // first wall tangetional unit vector
  CFreal lx = ny*vz - nz*vy;
  CFreal ly = nz*vx - nx*vz;
  CFreal lz = nx*vy - ny*vx;
  const CFreal labs = std::sqrt(lx*lx + ly*ly + lz*lz);

  lx = lx / labs;
  ly = ly / labs;
  lz = lz / labs;

  // second wall tangentional unit vector
  const CFreal mx = ny*lz - nz*ly;
  const CFreal my = nz*lx - nx*lz;
  const CFreal mz = nx*ly - ny*lx;
  
  const CFreal U = u*nx + v*ny + w*nz;
  const CFreal V = u*lx + v*ly + w*lz;
  const CFreal W = u*mx + v*my + w*mz;

  //cout << "n2  = " << nx*nx + ny*ny + nz*nz << endl;
  //cout << "v2  = " << vx*vx + vy*vy + vz*vz << endl;
  //cout << "l2  = " << lx*lx + ly*ly + lz*lz << endl;
  //cout << "m2  = " << mx*mx + my*my + mz*mz << endl;

  const CFreal H = lData[EulerTerm::H];
  const CFuint firstTv = getModel()->getFirstScalarVar(1);
  const CFreal ev = lData[firstTv];
  //  const CFreal a2 = (1. + beta)*p/rho;
  // const CFreal a = sqrt(a2);
  const CFreal a = lData[EulerTerm::A];
  const CFreal a2 = lData[EulerTerm::A]*lData[EulerTerm::A];

  const CFuint nbSpPlus1 = nbSpecies+1;
  const CFuint nbSpPlus2 = nbSpecies+2;
  const CFuint nbSpPlus3 = nbSpecies+3;
  const CFuint nbSpPlus4 = nbSpecies+4;

  // rightEv = 0.0;
  // leftEv = 0.0;
  for (CFuint is = 0; is < nbSpecies; ++is) {
    rightEv(is, is) = 1.;
    rightEv(is,nbSpPlus2) = 0.5*_ys[is];
    rightEv(is,nbSpPlus3) = 0.5*_ys[is];

    rightEv(nbSpecies,is) = u;
    rightEv(nbSpPlus1,is) = v;
    rightEv(nbSpPlus2,is) = w;
    rightEv(nbSpPlus3,is) = q - _alpha[is]/beta;
  }

  rightEv(nbSpecies,nbSpecies) = lx*a2;
  rightEv(nbSpecies,nbSpPlus1) = mx*a2;
  rightEv(nbSpecies,nbSpPlus2) = 0.5*(u+a*nx);
  rightEv(nbSpecies,nbSpPlus3) = 0.5*(u-a*nx);

  rightEv(nbSpPlus1,nbSpecies) = ly*a2;
  rightEv(nbSpPlus1,nbSpPlus1) = my*a2;
  rightEv(nbSpPlus1,nbSpPlus2) = 0.5*(v+a*ny); // 0.5(v+a*ny) ???
  rightEv(nbSpPlus1,nbSpPlus3) = 0.5*(v-a*ny); // 0.5(v+a*ny) ???

  rightEv(nbSpPlus2,nbSpecies) = lz*a2;
  rightEv(nbSpPlus2,nbSpPlus1) = mz*a2;
  rightEv(nbSpPlus2,nbSpPlus2) = 0.5*(w+a*nz);
  rightEv(nbSpPlus2,nbSpPlus3) = 0.5*(w-a*nz);

  rightEv(nbSpPlus3,nbSpecies) = V*a2;
  rightEv(nbSpPlus3,nbSpPlus1) = W*a2;
  rightEv(nbSpPlus3,nbSpPlus2) = 0.5*(H + a*U);
  rightEv(nbSpPlus3,nbSpPlus3) = 0.5*(H - a*U);
  rightEv(nbSpPlus3,nbSpPlus4) = -phi/beta;

  rightEv(nbSpPlus4,nbSpPlus2) = 0.5*ev;
  rightEv(nbSpPlus4,nbSpPlus3) = 0.5*ev;
  rightEv(nbSpPlus4,nbSpPlus4) = 1.0;

  rightEv /= a2;

  //  cout << "R = "<< endl;
  //   cout << rightEv << endl <<endl;

  const CFreal bq = beta*q;
  for (CFuint is = 0; is < nbSpecies; ++is) {
    for (CFuint js = 0; js < nbSpecies; ++js) {
      const CFreal a2delta = (js != is) ? 0.0 : a2;
      leftEv(is,js) = a2delta - _ys[is]*(_alpha[js] + bq);
    }

    leftEv(is,nbSpecies) = beta*u*_ys[is];
    leftEv(is,nbSpPlus1) = beta*v*_ys[is];
    leftEv(is,nbSpPlus2) = beta*w*_ys[is];
    leftEv(is,nbSpPlus3) = -beta*_ys[is];
    leftEv(is,nbSpPlus4) = beta*_ys[is];

    const CFreal alphaPlusBetaq = _alpha[is] + bq;
    leftEv(nbSpecies,is) = -V;
    leftEv(nbSpPlus1,is) = -W;
    leftEv(nbSpPlus2,is) = alphaPlusBetaq - U*a;
    leftEv(nbSpPlus3,is) = alphaPlusBetaq + U*a;
    leftEv(nbSpPlus4,is) = -ev*alphaPlusBetaq;
  }

  leftEv(nbSpecies,nbSpecies) = lx;
  leftEv(nbSpecies,nbSpPlus1) = ly;
  leftEv(nbSpecies,nbSpPlus2) = lz;

  leftEv(nbSpPlus1,nbSpecies) = mx;
  leftEv(nbSpPlus1,nbSpPlus1) = my;
  leftEv(nbSpPlus1,nbSpPlus2) = mz;

  leftEv(nbSpPlus2,nbSpecies) = a*nx - beta*u; //a*nx - beta*u;
  leftEv(nbSpPlus2,nbSpPlus1) = a*ny - beta*v; //a*ny - beta*v;
  leftEv(nbSpPlus2,nbSpPlus2) = a*nz - beta*w;
  leftEv(nbSpPlus2,nbSpPlus3) = beta;
  leftEv(nbSpPlus2,nbSpPlus4) = phi;

  leftEv(nbSpPlus3,nbSpecies) = -a*nx - beta*u;
  leftEv(nbSpPlus3,nbSpPlus1) = -a*ny - beta*v;
  leftEv(nbSpPlus3,nbSpPlus2) = -a*nz - beta*w;
  leftEv(nbSpPlus3,nbSpPlus3) = beta;
  leftEv(nbSpPlus3,nbSpPlus4) = phi;

  leftEv(nbSpPlus4,nbSpecies) = beta*u*ev;
  leftEv(nbSpPlus4,nbSpPlus1) = beta*v*ev;
  leftEv(nbSpPlus4,nbSpPlus2) = beta*w*ev;
  leftEv(nbSpPlus4,nbSpPlus3) = -beta*ev;
  leftEv(nbSpPlus4,nbSpPlus4) = a2 - phi*ev;

  eValues = U;
  eValues[nbSpPlus2] += a;
  eValues[nbSpPlus3] -= a;

  //  // DEBUG!!!

  // cout << "L = "<< endl;
  //   cout << _leftEv << endl <<endl;
  
  //  RealMatrix mat(rightEv.nbRows(), leftEv.nbRows());
  //   mat = _rightEv*_leftEv;
  //   cout <<"R*L" << endl;
  //   cout << mat << endl << endl;
  
  //   mat = _leftEv*_rightEv;
  //   cout <<"L*R" << endl;
  //   cout << mat << endl << endl;

  // jacobian matrix
  //  mat = 0.0;
  //  for (CFuint is = 0; is < nbSpecies; ++is) {
  //    for (CFuint js = 0; js < nbSpecies; ++js) {
  //      const CFreal delta = (js != is) ? 0.0 : 1.0; 
  //      mat(is,js) = (delta - _ys[is])*U;
  //    }

  //    mat(is,nbSpecies) = _ys[is]*nx;
  //    mat(is,nbSpPlus1) = _ys[is]*ny;
  //    mat(is,nbSpPlus2) = _ys[is]*nz;
  //    mat(is,nbSpPlus3) = 0.0;
  //    mat(is,nbSpPlus4) = 0.0;
    
  //    const CFreal alphaPlusBeta = _alpha[is] + bq;
  //    mat(nbSpecies,is) = alphaPlusBeta*nx -U*u;
  //    mat(nbSpPlus1,is) = alphaPlusBeta*ny -U*v;
  //    mat(nbSpPlus2,is) = alphaPlusBeta*nz -U*w;
  //    mat(nbSpPlus3,is) = (alphaPlusBeta - H)*U;
  //    mat(nbSpPlus4,is) = -U*ev;
  //  }

  //  mat(nbSpecies,nbSpecies) = (1.-beta)*u*nx + U;
  //  mat(nbSpecies,nbSpPlus1) = -beta*v*nx + u*ny;
  //  mat(nbSpecies,nbSpPlus2) = -beta*w*nx + u*nz;
  //  mat(nbSpecies,nbSpPlus3) = beta*nx;
  //  mat(nbSpecies,nbSpPlus4) = -beta*nx;
  
  //  mat(nbSpPlus1,nbSpecies) = -beta*u*ny + v*nx;
  //  mat(nbSpPlus1,nbSpPlus1) = (1.-beta)*v*ny + U;
  //  mat(nbSpPlus1,nbSpPlus2) = -beta*w*ny + v*nz;
  //  mat(nbSpPlus1,nbSpPlus3) = beta*ny;
  //  mat(nbSpPlus1,nbSpPlus4) = -beta*ny;
  
  //  mat(nbSpPlus2,nbSpecies) = -beta*u*nz + w*nx;
  //  mat(nbSpPlus2,nbSpPlus1) = -beta*v*nz + w*ny;
  //  mat(nbSpPlus2,nbSpPlus2) = (1.-beta)*w*nz + U;
  //  mat(nbSpPlus2,nbSpPlus3) = beta*nz;
  //  mat(nbSpPlus2,nbSpPlus4) = -beta*nz;

  //  mat(nbSpPlus3,nbSpecies) = H*nx - beta*U*u;
  //  mat(nbSpPlus3,nbSpPlus1) = H*ny - beta*U*v;
  //  mat(nbSpPlus3,nbSpPlus2) = H*nz - beta*U*w;
  //  mat(nbSpPlus3,nbSpPlus3) = (1.+beta)*U;
  //  mat(nbSpPlus3,nbSpPlus4) = -beta*U;
  
  //  mat(nbSpPlus4,nbSpecies) = ev*nx;
  //  mat(nbSpPlus4,nbSpPlus1) = ev*ny;
  //  mat(nbSpPlus4,nbSpPlus2) = ev*nz;
  //  mat(nbSpPlus4,nbSpPlus3) = 0.0;
  //  mat(nbSpPlus4,nbSpPlus4) = U;  
  
  //  RealMatrix m1(10,10);
  //  RealMatrix m2(10,10);
  
  //  m1 = mat*rightEv;
  //  m2 = leftEv*m1;
  
  //  cout << "lambda M = " << endl;
  //  cout << m2 << endl;
  //  cout << "lambda = " << eValues << endl << endl;

  //  for (CFuint i = 0; i < 10; ++i) {
  //    m2(i,i) = m2(i,i) - eValues[i];
  //  }

  //  CFreal matMax = 0.0;
  //  CFreal displayLimit = 1.0e-5;

  //  for (CFuint i = 0; i < 10; ++i) {
  //    for (CFuint j = 0; j < 10; ++j) {
  //      if ( (m2(i,j)*m2(i,j)) > (matMax*matMax) ) {
  //        matMax = m2(i,j);
  //      }
  //    }
  //  }

  //  if ( (matMax*matMax) >  (displayLimit*displayLimit) ) {
  //    cout << "Max difference between eigenvalues:  " << matMax << endl;
  //  }

  //  // DEBUG end!!!

}

//////////////////////////////////////////////////////////////////////////////

CFuint Euler3DNEQCons::getBlockSeparator() const
{
  return Framework::PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQCons::splitJacobian(RealMatrix& jacobPlus,
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
  const CFreal beta = p/(rho*T)/cvTr;
  const CFuint firstSpecies = getModel()->getFirstScalarVar(0);

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
  }
  
  CFreal phi = 0.0; 
  if (_library->presenceElectron()) {
    // assume that electrons have 0 as mixture ID  
    phi = _RiGas[0]*_ys[0]/eData->dEvTv - beta;
  }
  else {
    phi = -beta;
  }
  
  CFreal Tq = 0.0;
  for (CFuint is = 0; is < nbSpecies; ++is) {
    if (!_library->presenceElectron()) {
      Tq = T;
    }
    else {
      if (is == 0) {
	Tq = 0.0; // get vibrational temperature
      }
    }
    
    _alpha[is] = _RiGas[is]*Tq - beta*(eData->energyTr)[is]; // change here to include phi
  }
  
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal nz = normal[ZZ];
  const CFreal u = lData[EulerTerm::VX];
  const CFreal v = lData[EulerTerm::VY];
  const CFreal w = lData[EulerTerm::VZ];

  const CFreal V2 = lData[EulerTerm::V]*lData[EulerTerm::V];
  const CFreal q = 0.5*V2;

  // defining wall tangentional components

  // arbitrary unit vector (not parallel to unit normal)
  const CFreal vx = nz;
  const CFreal vy = nx;
  const CFreal vz = ny;

  // first wall tangetional unit vector
  CFreal lx = ny*vz - nz*vy;
  CFreal ly = nz*vx - nx*vz;
  CFreal lz = nx*vy - ny*vx;
  const CFreal labs = std::sqrt(lx*lx + ly*ly + lz*lz);

  lx = lx / labs;
  ly = ly / labs;
  lz = lz / labs;

  // second wall tangentional unit vector
  const CFreal mx = ny*lz - nz*ly;
  const CFreal my = nz*lx - nx*lz;
  const CFreal mz = nx*ly - ny*lx;
  
  const CFreal U = u*nx + v*ny + w*nz;
  const CFreal V = u*lx + v*ly + w*lz;
  const CFreal W = u*mx + v*my + w*mz;

  //cout << "n2  = " << nx*nx + ny*ny + nz*nz << endl;
  //cout << "v2  = " << vx*vx + vy*vy + vz*vz << endl;
  //cout << "l2  = " << lx*lx + ly*ly + lz*lz << endl;
  //cout << "m2  = " << mx*mx + my*my + mz*mz << endl;

  const CFreal H = lData[EulerTerm::H];
  const CFuint firstTv = getModel()->getFirstScalarVar(1);
  const CFreal ev = lData[firstTv];
  const CFreal a2 = (1. + beta)*p/rho;
  const CFreal a = sqrt(a2);

  const CFuint nbSpPlus1 = nbSpecies+1; 
  const CFuint nbSpPlus2 = nbSpecies+2; 
  const CFuint nbSpPlus3 = nbSpecies+3;
  const CFuint nbSpPlus4 = nbSpecies+4;

  // _rightEv = 0.0;
  // _leftEv = 0.0;
  for (CFuint is = 0; is < nbSpecies; ++is) {
    _rightEv(is, is) = 1.;
    _rightEv(is,nbSpPlus2) = 0.5*_ys[is];
    _rightEv(is,nbSpPlus3) = 0.5*_ys[is];

    _rightEv(nbSpecies,is) = u;
    _rightEv(nbSpPlus1,is) = v;
    _rightEv(nbSpPlus2,is) = w;
    _rightEv(nbSpPlus3,is) = q - _alpha[is]/beta;
  }  

  _rightEv(nbSpecies,nbSpecies) = lx*a2;
  _rightEv(nbSpecies,nbSpPlus1) = mx*a2;
  _rightEv(nbSpecies,nbSpPlus2) = 0.5*(u+a*nx);
  _rightEv(nbSpecies,nbSpPlus3) = 0.5*(u-a*nx);

  _rightEv(nbSpPlus1,nbSpecies) = ly*a2;
  _rightEv(nbSpPlus1,nbSpPlus1) = my*a2;
  _rightEv(nbSpPlus1,nbSpPlus2) = 0.5*(v+a*ny); // 0.5(v+a*ny) ???
  _rightEv(nbSpPlus1,nbSpPlus3) = 0.5*(v-a*ny); // 0.5(v+a*ny) ???

  _rightEv(nbSpPlus2,nbSpecies) = lz*a2;
  _rightEv(nbSpPlus2,nbSpPlus1) = mz*a2;
  _rightEv(nbSpPlus2,nbSpPlus2) = 0.5*(w+a*nz);
  _rightEv(nbSpPlus2,nbSpPlus3) = 0.5*(w-a*nz);

  _rightEv(nbSpPlus3,nbSpecies) = V*a2;
  _rightEv(nbSpPlus3,nbSpPlus1) = W*a2;
  _rightEv(nbSpPlus3,nbSpPlus2) = 0.5*(H + a*U);
  _rightEv(nbSpPlus3,nbSpPlus3) = 0.5*(H - a*U);
  _rightEv(nbSpPlus3,nbSpPlus4) = -phi/beta;

  _rightEv(nbSpPlus4,nbSpPlus2) = 0.5*ev;
  _rightEv(nbSpPlus4,nbSpPlus3) = 0.5*ev;
  _rightEv(nbSpPlus4,nbSpPlus4) = 1.0;
  
  _rightEv /= a2;

 //  cout << "R = "<< endl;
//   cout << _rightEv << endl <<endl;

  const CFreal bq = beta*q;
  for (CFuint is = 0; is < nbSpecies; ++is) {
    for (CFuint js = 0; js < nbSpecies; ++js) {
      const CFreal a2delta = (js != is) ? 0.0 : a2;
      _leftEv(is,js) = a2delta - _ys[is]*(_alpha[js] + bq);
    }

    _leftEv(is,nbSpecies) = beta*u*_ys[is];
    _leftEv(is,nbSpPlus1) = beta*v*_ys[is];
    _leftEv(is,nbSpPlus2) = beta*w*_ys[is];
    _leftEv(is,nbSpPlus3) = -beta*_ys[is];
    _leftEv(is,nbSpPlus4) = -phi*_ys[is];

    const CFreal alphaPlusBetaq = _alpha[is] + bq;
    _leftEv(nbSpecies,is) = -V;
    _leftEv(nbSpPlus1,is) = -W;
    _leftEv(nbSpPlus2,is) = alphaPlusBetaq - U*a;
    _leftEv(nbSpPlus3,is) = alphaPlusBetaq + U*a;
    _leftEv(nbSpPlus4,is) = -ev*alphaPlusBetaq;

  }  

  _leftEv(nbSpecies,nbSpecies) = lx;
  _leftEv(nbSpecies,nbSpPlus1) = ly;
  _leftEv(nbSpecies,nbSpPlus2) = lz;

  _leftEv(nbSpPlus1,nbSpecies) = mx;
  _leftEv(nbSpPlus1,nbSpPlus1) = my;
  _leftEv(nbSpPlus1,nbSpPlus2) = mz;

  _leftEv(nbSpPlus2,nbSpecies) = a*nx - beta*u; //a*nx - beta*u;
  _leftEv(nbSpPlus2,nbSpPlus1) = a*ny - beta*v; //a*ny - beta*v;
  _leftEv(nbSpPlus2,nbSpPlus2) = a*nz - beta*w;
  _leftEv(nbSpPlus2,nbSpPlus3) = beta;
  _leftEv(nbSpPlus2,nbSpPlus4) = phi;

  _leftEv(nbSpPlus3,nbSpecies) = -a*nx - beta*u;
  _leftEv(nbSpPlus3,nbSpPlus1) = -a*ny - beta*v;
  _leftEv(nbSpPlus3,nbSpPlus2) = -a*nz - beta*w;
  _leftEv(nbSpPlus3,nbSpPlus3) = beta;
  _leftEv(nbSpPlus3,nbSpPlus4) = phi;

  _leftEv(nbSpPlus4,nbSpecies) = beta*u*ev;
  _leftEv(nbSpPlus4,nbSpPlus1) = beta*v*ev;
  _leftEv(nbSpPlus4,nbSpPlus2) = beta*w*ev;
  _leftEv(nbSpPlus4,nbSpPlus3) = -beta*ev;
  _leftEv(nbSpPlus4,nbSpPlus4) = a2 - phi*ev;
  
  // cout << "L = "<< endl;
  //   cout << _leftEv << endl <<endl;
  
  //  RealMatrix mat(_rightEv.nbRows(), _leftEv.nbRows());
  //   mat = _rightEv*_leftEv;
  //   cout <<"R*L" << endl;
  //   cout << mat << endl << endl;
  
  //   mat = _leftEv*_rightEv;
  //   cout <<"L*R" << endl;
  //   cout << mat << endl << endl;
  
  
  eValues = U;
  eValues[nbSpPlus2] += a;
  eValues[nbSpPlus3] -= a;
  
  // jacobian matrix DEBUG
  //   mat = 0.0;
  //   for (CFuint is = 0; is < nbSpecies; ++is) {
  //     for (CFuint js = 0; js < nbSpecies; ++js) {
  //       const CFreal delta = (js != is) ? 0.0 : 1.0; 
  //       mat(is,js) = (delta - _ys[is])*U;
  //     }

  //     mat(is,nbSpecies) = _ys[is]*nx;
  //     mat(is,nbSpPlus1) = _ys[is]*ny;
  //     mat(is,nbSpPlus2) = _ys[is]*nz;
  //     mat(is,nbSpPlus3) = 0.0;
  //     mat(is,nbSpPlus4) = 0.0;
    
  //     const CFreal alphaPlusBeta = _alpha[is] + bq;
  //     mat(nbSpecies,is) = alphaPlusBeta*nx -U*u;
  //     mat(nbSpPlus1,is) = alphaPlusBeta*ny -U*v;
  //     mat(nbSpPlus2,is) = alphaPlusBeta*nz -U*w;
  //     mat(nbSpPlus3,is) = (alphaPlusBeta - H)*U;
  //     mat(nbSpPlus4,is) = -U*ev;
  //   }

  //   mat(nbSpecies,nbSpecies) = (1.-beta)*u*nx + U;
  //   mat(nbSpecies,nbSpPlus1) = -beta*v*nx + u*ny;
  //   mat(nbSpecies,nbSpPlus2) = -beta*w*nx + u*nz;
  //   mat(nbSpecies,nbSpPlus3) = beta*nx;
  //   mat(nbSpecies,nbSpPlus4) = -beta*nx;
  
  //   mat(nbSpPlus1,nbSpecies) = -beta*u*ny + v*nx;
  //   mat(nbSpPlus1,nbSpPlus1) = (1.-beta)*v*ny + U;
  //   mat(nbSpPlus1,nbSpPlus2) = -beta*w*ny + v*nz;
  //   mat(nbSpPlus1,nbSpPlus3) = beta*ny;
  //   mat(nbSpPlus1,nbSpPlus4) = -beta*ny;
  
  //   mat(nbSpPlus2,nbSpecies) = -beta*u*nz + w*nx;
  //   mat(nbSpPlus2,nbSpPlus1) = -beta*v*nz + w*ny;
  //   mat(nbSpPlus2,nbSpPlus2) = (1.-beta)*w*nz + U;
  //   mat(nbSpPlus2,nbSpPlus3) = beta*nz;
  //   mat(nbSpPlus2,nbSpPlus4) = -beta*nz;

  //   mat(nbSpPlus3,nbSpecies) = H*nx - beta*U*u;
  //   mat(nbSpPlus3,nbSpPlus1) = H*ny - beta*U*v;
  //   mat(nbSpPlus3,nbSpPlus2) = H*nz - beta*U*w;
  //   mat(nbSpPlus3,nbSpPlus3) = (1.+beta)*U;
  //   mat(nbSpPlus3,nbSpPlus4) = -beta*U;
  
  //   mat(nbSpPlus4,nbSpecies) = ev*nx;
  //   mat(nbSpPlus4,nbSpPlus1) = ev*ny;
  //   mat(nbSpPlus4,nbSpPlus2) = ev*nz;
  //   mat(nbSpPlus4,nbSpPlus3) = 0.0;
  //   mat(nbSpPlus4,nbSpPlus4) = U;  
  
  //   RealMatrix m1(10,10);
  //   RealMatrix m2(10,10);
  
  //   m1 = mat*_rightEv;
  //   m2 = _leftEv*m1;
  
  //   cout << "lambda M = " << endl;
  //   cout << m2 << endl;
  //   cout << "lambda = " << eValues << endl << endl;
  
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
  CFLogDebugMax( "RightEigenvectors @Euler3DNEQCons::splitJacobian" << "\n"
		 << _rightEv << "\n");
  CFLogDebugMax( "LeftEigenvectors @Euler3DNEQCons::splitJacobian" << "\n"
		 << _leftEv << "\n");
  CFLogDebugMax( "EigenValues @Euler3DNEQCons::splitJacobian" << "\n"
		 << eValues << "\n" << "\n"); 
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQCons::computePhysicalData(const State& state,
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

  // U, V  and W Velocity Average
  data[EulerTerm::VX] = state[nbSpecies]*ovRho;
  data[EulerTerm::VY] = state[nbSpecies+1]*ovRho;
  data[EulerTerm::VZ] = state[nbSpecies+2]*ovRho;
  
  const CFreal V2 = data[EulerTerm::VX]*data[EulerTerm::VX] +
                    data[EulerTerm::VY]*data[EulerTerm::VY] +
                    data[EulerTerm::VZ]*data[EulerTerm::VZ];
  data[EulerTerm::V] = sqrt(V2);
  
  const CFuint startTv = getModel()->getFirstScalarVar(1);
  data[startTv]= state[nbSpecies+4]*ovRho; //ev

  data[EulerTerm::E]  = state[nbSpecies+3]*ovRho;
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

void Euler3DNEQCons::computeStateFromPhysicalData(const RealVector& data,
						 State& state)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DNEQCons::computeStateFromPhysicalData()");

  //  state[0] = data[EulerTerm::RHO];
  //  state[1] = data[EulerTerm::RHO]*data[EulerTerm::VX];
  //  state[2] = data[EulerTerm::RHO]*data[EulerTerm::VY];
  //  state[3] = data[EulerTerm::RHO]*data[EulerTerm::VZ];
  //  state[4] = data[EulerTerm::RHO]*data[EulerTerm::E];
  
  //  // Set the species
  //  const CFuint firstSpecies = getModel()->getFirstScalarVar(0);
  //  const CFuint nbSpecies = getModel()->getNbScalarVars(0);

  //  for (CFuint ie = 0; ie < nbSpecies; ++ie){
  //    state[5 + ie] = data[EulerTerm::RHO]*data[firstSpecies + ie];
  //  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal Euler3DNEQCons::getSpeed(const State& state) const
{
  const CFreal u = state[1]/state[0];
  const CFreal v = state[2]/state[0];
  const CFreal w = state[3]/state[0];
  return sqrt(u*u + v*v + w*w);
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQCons::setDimensionalValues(const State& state,
                                          RealVector& result)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DNEQCons::setDimensionalValues()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQCons::setAdimensionalValues(const State& state,
                                          RealVector& result)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DNEQCons::setAdimensionalValues()");
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQCons::setDimensionalValuesPlusExtraValues
(const State& state, RealVector& result,
 RealVector& extra)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DNEQCons::setDimensionalValuesPlusExtraValues()");
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> Euler3DNEQCons::getExtraVarNames() const
{
  cf_assert(_library.isNotNull());

  vector<std::string> names(3);
  names[0] = "rho";
  names[1] = "H";
  names[2] = "M";

  return names;
}

//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQCons::setup()
{
  MultiScalarVarSet<Euler3DVarSet>::setup();
  
  const CFuint nbSpecies = getModel()->getNbScalarVars(0);
  // set the equation set data for each of the equation subsets
  // first equation subset
  MultiScalarVarSet<Euler3DVarSet>::getEqSetData().resize(2);
  MultiScalarVarSet<Euler3DVarSet>::getEqSetData()[0].setup(0,0,nbSpecies);
  
  // second equation subset
  Euler3DVarSet::getEqSetData().resize(1);
  Euler3DVarSet::getEqSetData()[0].setup(1,nbSpecies,4);
  
  const CFuint nbTv = getModel()->getNbScalarVars(1);
  // third equation subset
  MultiScalarVarSet<Euler3DVarSet>::getEqSetData()[1].setup(2,nbSpecies + 4,nbTv);
  
  _library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(_library.isNotNull());
  
  _Rgas = _library->getRgas();
  
  _dhe.resize(3 + nbTv);
  _ys.resize(nbSpecies);
  const CFuint totNbEqs = 4 + nbTv + nbSpecies;
   
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

void Euler3DNEQCons::computePerturbedPhysicalData(const Framework::State& state,
						  const RealVector& pdataBkp,
						  RealVector& pdata,
						  CFuint iVar)
{
  throw Common::NotImplementedException (FromHere(),"Euler3DNEQCons::computePerturbedPhysicalData()");
}
      
//////////////////////////////////////////////////////////////////////////////

void Euler3DNEQCons::computeProjectedJacobian(const RealVector& normal,
				 RealMatrix& jacob)
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
  const CFreal nz = normal[ZZ];
  const CFreal u = lData[EulerTerm::VX];
  const CFreal v = lData[EulerTerm::VY];
  const CFreal w = lData[EulerTerm::VZ];

  const CFreal V2 = lData[EulerTerm::V]*lData[EulerTerm::V];
  const CFreal q = 0.5*V2;

  // defining wall tangentional components

  // arbitrary unit vector (not parallel to unit normal)
  const CFreal vx = nz;
  const CFreal vy = nx;
  const CFreal vz = ny;

  // first wall tangetional unit vector
  CFreal lx = ny*vz - nz*vy;
  CFreal ly = nz*vx - nx*vz;
  CFreal lz = nx*vy - ny*vx;
  const CFreal labs = std::sqrt(lx*lx + ly*ly + lz*lz);

  lx = lx / labs;
  ly = ly / labs;
  lz = lz / labs;

  // second wall tangentional unit vector
  // const CFreal mx = ny*lz - nz*ly;
  // const CFreal my = nz*lx - nx*lz;
  // const CFreal mz = nx*ly - ny*lx;
  
  const CFreal U = u*nx + v*ny + w*nz;
  // const CFreal V = u*lx + v*ly + w*lz;
  // const CFreal W = u*mx + v*my + w*mz;

  //cout << "n2  = " << nx*nx + ny*ny + nz*nz << endl;
  //cout << "v2  = " << vx*vx + vy*vy + vz*vz << endl;
  //cout << "l2  = " << lx*lx + ly*ly + lz*lz << endl;
  //cout << "m2  = " << mx*mx + my*my + mz*mz << endl;

  const CFreal H = lData[EulerTerm::H];
  const CFuint firstTv = getModel()->getFirstScalarVar(1);
  const CFreal ev = lData[firstTv];
  // const CFreal a2 = (1. + beta)*p/rho;
  const CFreal bq = beta*q;

  const CFuint nbSpPlus1 = nbSpecies+1; 
  const CFuint nbSpPlus2 = nbSpecies+2; 
  const CFuint nbSpPlus3 = nbSpecies+3;
  const CFuint nbSpPlus4 = nbSpecies+4;
  
  for (CFuint is = 0; is < nbSpecies; ++is) {
    for (CFuint js = 0; js < nbSpecies; ++js) {
      const CFreal delta = (js != is) ? 0.0 : 1.0; 
      jacob(is,js) = (delta - _ys[is])*U;
    }
    jacob(is,nbSpecies) = _ys[is]*nx;
    jacob(is,nbSpPlus1) = _ys[is]*ny;
    jacob(is,nbSpPlus2) = _ys[is]*nz;
    jacob(is,nbSpPlus3) = 0.0;
    jacob(is,nbSpPlus4) = 0.0;
    
    const CFreal alphaPlusBeta = _alpha[is] + bq;
    jacob(nbSpecies,is) = alphaPlusBeta*nx -U*u;
    jacob(nbSpPlus1,is) = alphaPlusBeta*ny -U*v;
    jacob(nbSpPlus2,is) = alphaPlusBeta*nz -U*w;
    jacob(nbSpPlus3,is) = (alphaPlusBeta - H)*U;
    jacob(nbSpPlus4,is) = -U*ev;
  }
  
  jacob(nbSpecies,nbSpecies) = (1.-beta)*u*nx + U;
  jacob(nbSpecies,nbSpPlus1) = -beta*v*nx + u*ny;
  jacob(nbSpecies,nbSpPlus2) = -beta*w*nx + u*nz;
  jacob(nbSpecies,nbSpPlus3) = beta*nx;
  jacob(nbSpecies,nbSpPlus4) = -beta*nx;
  
  jacob(nbSpPlus1,nbSpecies) = -beta*u*ny + v*nx;
  jacob(nbSpPlus1,nbSpPlus1) = (1.-beta)*v*ny + U;
  jacob(nbSpPlus1,nbSpPlus2) = -beta*w*ny + v*nz;
  jacob(nbSpPlus1,nbSpPlus3) = beta*ny;
  jacob(nbSpPlus1,nbSpPlus4) = -beta*ny;
  
  jacob(nbSpPlus2,nbSpecies) = -beta*u*nz + w*nx;
  jacob(nbSpPlus2,nbSpPlus1) = -beta*v*nz + w*ny;
  jacob(nbSpPlus2,nbSpPlus2) = (1.-beta)*w*nz + U;
  jacob(nbSpPlus2,nbSpPlus3) = beta*nz;
  jacob(nbSpPlus2,nbSpPlus4) = -beta*nz;

  jacob(nbSpPlus3,nbSpecies) = H*nx - beta*U*u;
  jacob(nbSpPlus3,nbSpPlus1) = H*ny - beta*U*v;
  jacob(nbSpPlus3,nbSpPlus2) = H*nz - beta*U*w;
  jacob(nbSpPlus3,nbSpPlus3) = (1.+beta)*U;
  jacob(nbSpPlus3,nbSpPlus4) = -beta*U;
  
  jacob(nbSpPlus4,nbSpecies) = ev*nx;
  jacob(nbSpPlus4,nbSpPlus1) = ev*ny;
  jacob(nbSpPlus4,nbSpPlus2) = ev*nz;
  jacob(nbSpPlus4,nbSpPlus3) = 0.0;
  jacob(nbSpPlus4,nbSpPlus4) = U;
}

//////////////////////////////////////////////////////////////////////////////
      
void Euler3DNEQCons::setStateVelocityIDs (std::vector<CFuint>& velIDs) 
{
  const CFuint nbSpecies = _library->getNbSpecies();
  velIDs.resize(3); 
  velIDs[XX] = nbSpecies; velIDs[YY] = nbSpecies + 1; velIDs[ZZ] = nbSpecies + 2;
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
