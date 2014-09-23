#include "Common/COOLFluiD.hh"

#include "MathTools/MathFunctions.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"

#include "DiscontGalerkin/StdBaseSolve.hh"
#include "DiscontGalerkin/DiscontGalerkin.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StdBaseSolve,DiscontGalerkinSolverData,DiscontGalerkinModule >
  stdBaseSolveProvider("StdBaseSolve");

//////////////////////////////////////////////////////////////////////////////

StdBaseSolve::StdBaseSolve(const std::string& name)
  : DiscontGalerkinSolverCom(name)
{
}

//////////////////////////////////////////////////////////////////////////////

StdBaseSolve::~StdBaseSolve()
{
}

//////////////////////////////////////////////////////////////////////////////

void StdBaseSolve::setup()
{
  CFAUTOTRACE;
  m_state = new Framework::State;
}

//////////////////////////////////////////////////////////////////////////////

void StdBaseSolve::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

// void StdBaseSolve::compute_EigenValVec_old_2D(Framework::State state, RealMatrix T, RealMatrix T1, RealMatrix *Pplus, RealMatrix *Pminus, RealMatrix *EigenVal, RealVector normal)
// {
// CFAUTOTRACE;
//   CFreal kap = 1.4;
//   CFreal kap1 = 0.4;
//   RealMatrix Q(4,4);
//   RealMatrix Q1(4,4);
//   RealMatrix TT(4,4);
//   RealMatrix TT1(4,4);
// //     !  r-density, u,v-velocity, p-pressure, kappa-poisson constant
// //     !  c-the local speed of sound, h-enthalpy
//
//   CFreal rlen = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
//
//   CFreal normalv[2] = {normal[0]/rlen,normal[1]/rlen};
//
//   //CFreal n12 = sqrt(normalv[0]*normalv[0] + normalv[1]*normalv[1]);
//
// // if ((normalv[0] != normal[0]/rlen) || (normalv[1] != normal[1]/rlen)) CFout << " ERROR " << CFendl;

//   CFreal r=state[0];
//   cf_assert(r > 0);
//
//   CFreal u = (state[1]*normalv[0] + state[2]*normalv[1])/r;
//   CFreal v = (-state[1]*normalv[1] + state[2]*normalv[0])/r;
//   CFreal vv=(u*u + v*v)/2.0;
//   CFreal p=kap1*(state[3]-r*vv);
//   cf_assert(p > 0);
//
//   CFreal c=sqrt(kap*p/r);
//   CFreal c2=c*c;
// //      T and T1 - transformation matrices
// for(CFuint i=0;i<4;i++)
//   for(CFuint j=0;j<4;j++)
//   {
//   Q(i,j)=0.0;
//   Q1(i,j)=0.0;
//   }
//
// Q(0,0) = 1.0;
// Q(1,1) = normalv[0];
// Q(1,2) = normalv[1];
// Q(2,1) = -normalv[1];
// Q(2,2) = normalv[0];
// Q(3,3) = 1.0;
// Q1(0,0) = 1.0;
// Q1(1,1) = 1.0;
// Q1(2,2) = 1.0;
// Q1(3,3) = 1.0;
//
//   for(CFuint i = 1; i<3;i++)
//   {
//     if (Q(i,i) == 0.0)
//     {
//       CFuint j=i+1;
//       while ((Q(j,i) == 0.0)&&(j<3)) j++;
//       if (j == 3)
//       {
//         cf_assert(Q(i,i)!=0.0);
//       }
//       else
//       {
//         for(CFuint k=1;k<3;k++)
//         {
//           swap(Q(i,k),Q(j,k));
//           swap(Q1(i,k),Q1(j,k));
//         }
//       }
//     }
//     for(CFuint j=2;j>0;j--)
//     {
//       Q1(i,j) /=Q(i,i);
//     }
//     for(CFuint j=2;j>i-1;j--)
//     {
//       Q(i,j) /=Q(i,i);
//     }
//     for(CFuint j = 1; j<3;j++)
//     {
//       if (j!=i)
//       {
//         for(CFuint k=2;k>0;k--)
//         {
//           Q1(j,k) -= Q1(i,k)*Q(j,i);
//         }
//         for(CFuint k=2;k>i-1;k--)
//         {
//           Q(j,k) -= Q(i,k)*Q(j,i);
//         }
//       }
//     }
//   }
//
// Q(0,0) = 1.0;
// Q(1,1) = normalv[0];
// Q(1,2) = normalv[1];
// Q(2,1) = -normalv[1];
// Q(2,2) = normalv[0];
// Q(3,3) = 1.0;
//
//    T(0,0) = 1.;
//    T(0,1) = 1.;
//    T(0,2) = 1.;
//    T(0,3) = 1.;
//    T(1,0) = u - c;
//    T(1,1) = u;
//    T(1,2) = u;
//    T(1,3) = u + c;
//    T(2,0) = v;
//    T(2,1) = v;
//    T(2,2) = v - c;
//    T(2,3) = v;
//    T(3,0) = vv + c2/kap1 - u*c;
//    T(3,1) = vv;
//    T(3,2) = vv - v*c;
//    T(3,3) = vv + c2/kap1 + u*c;
//
//     T1(0,0) = (kap1*vv + u*c)/2.0;
//     T1(0,1) = -(c + kap1*u)/2.0;
//     T1(0,2) = -v*kap1/2.0;
//     T1(0,3) = kap1/2.0;
//     T1(1,0) = c2 - c*v - kap1*vv;
//     T1(1,1) = u*kap1;
//     T1(1,2) = c + v*kap1;
//     T1(1,3) = -kap1;
//     T1(2,0) = v*c;
//     T1(2,1) = 0.0;
//     T1(2,2) = -c;
//     T1(2,3) = 0.0;
//     T1(3,0) = (kap1*vv - u*c)/2.0;
//     T1(3,1) = (c - kap1*u)/2.0;
//     T1(3,2) = -v*kap1/2.0;
//     T1(3,3) = kap1/2.0;
//
//     T1 /=c2;
//
// for(CFuint i=0;i<4;i++)
//   for(CFuint j=0;j<4;j++)
//   {
//   TT(i,j)=0.0;
//   TT1(i,j)=0.0;
//   }
//
// for(CFuint i=0;i<4;i++)
//   for(CFuint j=0;j<4;j++)
//     for(CFuint k=0;k<4;k++)
//     {
//       TT(i,j)+=Q1(i,k)*T(k,j);
//       TT1(i,j)+=T1(i,k)*Q(k,j);
//     }
//
//     // EigenVal ... eigenvalues of the matrix n1*a[w]+n2*b[w]
//
//     (*EigenVal)(0,1)=u;
//     (*EigenVal)(0,2)=(*EigenVal)(0,1);
//     (*EigenVal)(0,0)=(*EigenVal)(0,1)-c;
//     (*EigenVal)(0,3)=(*EigenVal)(0,1)+c;
//
//   for(CFuint i=0;i<4;i++)
//   {
//     if ((*EigenVal)(0,i) > 0)
//     {
//       (*EigenVal)(1,i)=0;
//     }
//     else
//     {
//       (*EigenVal)(1,i)=(*EigenVal)(0,i);
//       (*EigenVal)(0,i)=0;
//     }
//   }
// for(CFuint row=0;row<4;row++)
//   for(CFuint col=0;col<4;col++)
//   {
//     (*Pplus)(row,col)=0.0;
//     (*Pminus)(row,col)=0.0;
//     for(CFuint i=0;i<4;i++)
//     {
//       (*Pplus)(row,col)+=(*EigenVal)(0,i)*TT(row,i)*TT1(i,col);
//       (*Pminus)(row,col)+=(*EigenVal)(1,i)*TT(row,i)*TT1(i,col);
//     }
//   }
//   return;
// }

//////////////////////////////////////////////////////////////////////////////

CFuint StdBaseSolve::compute_EigenValVec2D(Framework::State state, RealMatrix T, RealMatrix T1, RealMatrix *Pplus, RealMatrix *Pminus, RealMatrix *EigenVal, RealVector normal)
{
CFAUTOTRACE;
  CFreal kap = 1.4;
  CFreal kap1 = 0.4;

//     !  r-density, u,v-velocity, p-pressure, kappa-poisson constant
//     !  c-the local speed of sound, h-enthalpy
  CFreal rlen = sqrt(normal[0]*normal[0] + normal[1]*normal[1] );
  CFreal normalv[2] = {normal[0]/rlen,normal[1]/rlen};

// if ((normalv[0] != normal[0]/rlen) || (normalv[1] != normal[1]/rlen)) CFout << " ERROR " << CFendl;

  CFreal r=state[0];
  cf_assert(r > 0);
  CFreal u=state[1]/r;
  CFreal v=state[2]/r;
  CFreal p=kap1*(state[3]-0.5*r*(u*u+v*v));
//   cf_assert(p > 0);
  if (p<=0)
  {
    return 1;
  }

     CFreal c=sqrt(kap*p/r);
     CFreal h=c*c/kap1+(u*u+v*v)/2.;
     CFreal c2=c*c;
     CFreal c22=2.*c2;
//      T and T1 - transformation matrices

   T(0,0) = 1.;
   T(0,1) = 0.;
   T(0,2) = 1./c22;
   T(0,3) = 1./c22;
   T(1,0) = u;
   T(1,1) = normalv[1];
   T(1,2) = (u+c*normalv[0])/c22;
   T(1,3) = (u-c*normalv[0])/c22;
   T(2,0) = v;
   T(2,1) = -normalv[0];
   T(2,2) = (v+c*normalv[1])/c22;
   T(2,3) = (v-c*normalv[1])/c22;
   T(3,0) = (u*u+v*v)/2.;
   T(3,1) = normalv[1]*u-normalv[0]*v;
   T(3,2) = (h+c*(normalv[0]*u+normalv[1]*v))/c22;
   T(3,3) = (h-c*(normalv[0]*u+normalv[1]*v))/c22;

    T1(0,0) = 1-kap1/c2*(u*u/2+v*v/2);
    T1(0,1) = kap1/c2*u;
    T1(0,2) = kap1/c2*v;
    T1(0,3) = -kap1/c2;
    T1(1,0) = normalv[0]*v-normalv[1]*u;
    T1(1,1) = normalv[1];
    T1(1,2) = -normalv[0];
    T1(1,3) = 0.;
    T1(2,0) = -c*(normalv[0]*u+normalv[1]*v)+kap1*(u*u/2+v*v/2);
    T1(2,1) = c*normalv[0]-kap1*u;
    T1(2,2) = c*normalv[1]-kap1*v;
    T1(2,3) = kap1;
    T1(3,0) = c*(normalv[0]*u+normalv[1]*v)+kap1*(u*u/2+v*v/2);
    T1(3,1) = -c*normalv[0]-kap1*u;
    T1(3,2) = -c*normalv[1]-kap1*v;
    T1(3,3) = kap1;


    // EigenVal ... eigenvalues of the matrix n1*a[w]+n2*b[w]

    (*EigenVal)(0,0)=normalv[0]*u+normalv[1]*v;
    (*EigenVal)(0,1)=(*EigenVal)(0,0);
    (*EigenVal)(0,2)=(*EigenVal)(0,0)+c;
    (*EigenVal)(0,3)=(*EigenVal)(0,0)-c;

// DataHandle<CFreal> maxEigen = socket_maxEigen.getDataHandle();
//     if ( maxEigen < EigenVal(0,2) )
//     {
//       maxEigen = EigenVal(0,2);
// CFout << "\n" << maxEigen << CFendl;
//     }

  for(CFuint i=0;i<4;i++)
  {
    if ((*EigenVal)(0,i) > 0)
    {
      (*EigenVal)(1,i)=0;
    }
    else
    {
      (*EigenVal)(1,i)=(*EigenVal)(0,i);
      (*EigenVal)(0,i)=0;
    }
  }
for(CFuint row=0;row<4;row++)
  for(CFuint col=0;col<4;col++)
  {
    (*Pplus)(row,col)=0.0;
    (*Pminus)(row,col)=0.0;
    for(CFuint i=0;i<4;i++)
    {
      (*Pplus)(row,col)+=(*EigenVal)(0,i)*T(row,i)*T1(i,col);
      (*Pminus)(row,col)+=(*EigenVal)(1,i)*T(row,i)*T1(i,col);
    }
  }
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

CFuint StdBaseSolve::compute_EigenValVec3D(Framework::State state, RealMatrix T, RealMatrix T1, RealMatrix *Pplus, RealMatrix *Pminus, RealMatrix *EigenVal, RealVector normal)
{
CFAUTOTRACE;
  CFreal kap = 1.4;
  CFreal kap1 = 0.4;
  RealMatrix Q(5,5);
  RealMatrix Q1(5,5);
  RealMatrix TT(5,5);
  RealMatrix TT1(5,5);

//     !  r-density, u,v-velocity, p-pressure, kappa-poisson constant
//     !  c-the local speed of sound, h-enthalpy
  CFreal rlen = sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);

  CFreal normalv[3] = {normal[0]/rlen,normal[1]/rlen,normal[2]/rlen};

  CFreal n12 = sqrt(normalv[0]*normalv[0] + normalv[1]*normalv[1]);

  CFreal r=state[0];
  cf_assert(r > 0);
   CFreal u;
   CFreal v;
   CFreal w;
  if (n12)
  {
    u = (state[1]*normalv[0] + state[2]*normalv[1] + state[3]*normalv[2])/r;
    v = (-state[1]*normalv[1] + state[2]*normalv[0])/(r*n12);
    w = (-state[1]*normalv[0]*normalv[2] - state[2]*normalv[1]*normalv[2] + state[3]*n12*n12)/(r*n12);
  }
  else
  {
    u = MathTools::MathFunctions::sign(normalv[2])*state[3]/r;
    v = state[2]/r;
    w = -MathTools::MathFunctions::sign(normalv[2])*state[1]/r;
  }
  CFreal vv=(u*u + v*v + w*w)/2.0;

  CFreal p=kap1*(state[4]-(state[1]*state[1]+state[2]*state[2]+state[3]*state[3])/(2.0*r));
//   cf_assert(p > 0);
  if (p<=0)
  {
    return 1;
  }

  CFreal c=sqrt(kap*p/r);
  CFreal c2=c*c;
//      T and T1 - transformation matrices
  for(CFuint i=0;i<5;i++)
    for(CFuint j=0;j<5;j++)
    {
      Q(i,j)=0.0;
      Q1(i,j)=0.0;
    }

  if (n12)
  {
    Q(0,0) = 1.0;
    Q(1,1) = normalv[0];
    Q(1,2) = normalv[1];
    Q(1,3) = normalv[2];
    Q(2,1) = -normalv[1]/n12;
    Q(2,2) = normalv[0]/n12;
    Q(3,1) = -normalv[0]*normalv[2]/n12;
    Q(3,2) = -normalv[1]*normalv[2]/n12;
    Q(3,3) = n12;
    Q(4,4) = 1.0;
  }
  else
  {
    Q(0,0) = 1.0;
    Q(1,1) = 0.0;
    Q(1,2) = 0.0;
    Q(1,3) = MathTools::MathFunctions::sign(normalv[2]);
    Q(2,1) = 0.0;
    Q(2,2) = 1.0;
    Q(3,1) = -MathTools::MathFunctions::sign(normalv[2]);
    Q(3,2) = 0.0;
    Q(3,3) = 0.0;
    Q(4,4) = 1.0;
  }

  Q1(0,0) = 1.0;
  Q1(1,1) = 1.0;
  Q1(2,2) = 1.0;
  Q1(3,3) = 1.0;
  Q1(4,4) = 1.0;

  for(CFuint i = 1; i<4;i++)
  {
    if (Q(i,i) == 0.0)
    {
      CFuint j=i+1;
      while ((Q(j,i) == 0.0)&&(j<4)) j++;
      if (j == 4)
      { cf_assert (Q(i,i)!=0.0); }
      else
      {
        for(CFuint k=1;k<4;k++)
        {
          swap(Q(i,k),Q(j,k));
          swap(Q1(i,k),Q1(j,k));
        }
      }
    }
    for(CFuint j=3;j>0;j--)
    {
      Q1(i,j) /=Q(i,i);
    }
    for(CFuint j=3;j>i-1;j--)
    {
      Q(i,j) /=Q(i,i);
    }
    for(CFuint j = 1; j<4;j++)
    {
      if (j!=i)
      {
        for(CFuint k=3;k>0;k--)
        {
          Q1(j,k) -= Q1(i,k)*Q(j,i);
        }
        for(CFuint k=3;k>i-1;k--)
        {
          Q(j,k) -= Q(i,k)*Q(j,i);
        }
      }
    }
  }


  if (n12)
  {
    Q(0,0) = 1.0;
    Q(1,1) = normalv[0];
    Q(1,2) = normalv[1];
    Q(1,3) = normalv[2];
    Q(2,1) = -normalv[1]/n12;
    Q(2,2) = normalv[0]/n12;
    Q(3,1) = -normalv[0]*normalv[2]/n12;
    Q(3,2) = -normalv[1]*normalv[2]/n12;
    Q(3,3) = n12;
    Q(4,4) = 1.0;
  }
// {
// Q(0,0) = 1.0;
// Q(1,1) = normalv[0];
// Q(1,2) = normalv[1];
// Q(1,3) = normalv[2];
// Q(2,1) = -MathTools::MathFunctions::sign(normalv[1])*normalv[1]/n12;
// Q(2,2) = MathTools::MathFunctions::sign(normalv[0])*normalv[0]/n12;
// Q(3,1) = -MathTools::MathFunctions::sign(normalv[0])*normalv[0]*normalv[2]/n12;
// Q(3,2) = -MathTools::MathFunctions::sign(normalv[1])*normalv[1]*normalv[2]/n12;
// Q(3,3) = n12;
// Q(4,4) = 1.0;
// }
  else
  {
    Q(0,0) = 1.0;
    Q(1,1) = 0.0;
    Q(1,2) = 0.0;
    Q(1,3) = MathTools::MathFunctions::sign(normalv[2]);
    Q(2,1) = 0.0;
    Q(2,2) = 1.0;
    Q(3,1) = -MathTools::MathFunctions::sign(normalv[2]);
    Q(3,2) = 0.0;
    Q(3,3) = 0.0;
    Q(4,4) = 1.0;
  }

  T(0,0) = 1.;
  T(0,1) = 1.;
  T(0,2) = 1.;
  T(0,3) = 1.;
  T(0,4) = 1.;
  T(1,0) = u - c;
  T(1,1) = u;
  T(1,2) = u;
  T(1,3) = u;
  T(1,4) = u + c;
  T(2,0) = v;
  T(2,1) = v;
  T(2,2) = v - c;
  T(2,3) = v;
  T(2,4) = v;
  T(3,0) = w;
  T(3,1) = w;
  T(3,2) = w;
  T(3,3) = w - c;
  T(3,4) = w;
  T(4,0) = vv + c2/kap1 - u*c;
  T(4,1) = vv;
  T(4,2) = vv - v*c;
  T(4,3) = vv - w*c;
  T(4,4) = vv + c2/kap1 + u*c;

  T1(0,0) = (kap1*vv + u*c)/2.0;
  T1(0,1) = -(c + kap1*u)/2.0;
  T1(0,2) = -v*kap1/2.0;
  T1(0,3) = -w*kap1/2.0;
  T1(0,4) = kap1/2.0;
  T1(1,0) = c2 - c*(v + w) - kap1*vv;
  T1(1,1) = u*kap1;
  T1(1,2) = c + v*kap1;
  T1(1,3) = c + w*kap1;
  T1(1,4) = -kap1;
  T1(2,0) = v*c;
  T1(2,1) = 0.0;
  T1(2,2) = -c;
  T1(2,3) = 0.0;
  T1(2,4) = 0.0;
  T1(3,0) = w*c;
  T1(3,1) = 0.0;
  T1(3,2) = 0.0;
  T1(3,3) = -c;
  T1(3,4) = 0.0;
  T1(4,0) = (kap1*vv - u*c)/2.0;
  T1(4,1) = (c - kap1*u)/2.0;
  T1(4,2) = -v*kap1/2.0;
  T1(4,3) = -w*kap1/2.0;
  T1(4,4) = kap1/2.0;

  T1 /=c2;

  for(CFuint i=0;i<5;i++)
    for(CFuint j=0;j<5;j++)
    {
      TT(i,j)=0.0;
      TT1(i,j)=0.0;
    }

  for(CFuint i=0;i<5;i++)
    for(CFuint j=0;j<5;j++)
      for(CFuint k=0;k<5;k++)
      {
        TT(i,j)+=Q1(i,k)*T(k,j);
        TT1(i,j)+=T1(i,k)*Q(k,j);
      }

    // EigenVal ... eigenvalues of the matrix n1*a[w]+n2*b[w]

  (*EigenVal)(0,1)=u;
  (*EigenVal)(0,2)=(*EigenVal)(0,1);
  (*EigenVal)(0,3)=(*EigenVal)(0,1);
  (*EigenVal)(0,0)=(*EigenVal)(0,1)-c;
  (*EigenVal)(0,4)=(*EigenVal)(0,1)+c;

  for(CFuint i=0;i<5;i++)
  {
    if ((*EigenVal)(0,i) > 0)
    {
      (*EigenVal)(1,i)=0;
    }
    else
    {
      (*EigenVal)(1,i)=(*EigenVal)(0,i);
      (*EigenVal)(0,i)=0;
    }
  }
  for(CFuint row=0;row<5;row++)
    for(CFuint col=0;col<5;col++)
    {
      (*Pplus)(row,col)=0.0;
      (*Pminus)(row,col)=0.0;
      for(CFuint i=0;i<5;i++)
      {
        (*Pplus)(row,col)+=(*EigenVal)(0,i)*TT(row,i)*TT1(i,col);
        (*Pminus)(row,col)+=(*EigenVal)(1,i)*TT(row,i)*TT1(i,col);
      }
    }
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

void StdBaseSolve::compute_Amatrix2D(Framework::State state, std::vector < RealMatrix > *m_A)
{
  CFreal gamma = 1.4;

  CFreal temp = (gamma - 1.0)*(state[1]*state[1] + state[2]*state[2]);

//   m_A[0][0][0]=0;
  (*m_A)[0](0,1)=1;
//   m_A[0][0][2]=0;
//   m_A[0][0][3]=0;

  (*m_A)[0](1,0)=-(state[1]*state[1] - 1.0/2.0*temp)/(state[0]*state[0]);
  (*m_A)[0](1,1)=(3.0 - gamma)*state[1]/state[0];
  (*m_A)[0](1,2)=(1.0 - gamma)*state[2]/state[0];
  (*m_A)[0](1,3)=gamma - 1;

  (*m_A)[0](2,0)=-state[1]*state[2]/(state[0]*state[0]);
  (*m_A)[0](2,1)=state[2]/state[0];
  (*m_A)[0](2,2)=state[1]/state[0];
//   (*m_A)[0][2][3]=0;

  (*m_A)[0](3,0)=state[1]/(state[0]*state[0])*(-gamma*state[3] + temp/state[0]);
  (*m_A)[0](3,1)=(gamma*state[3] - (temp/2.0 + (gamma - 1)*state[1]*state[1])/state[0])/state[0];
  (*m_A)[0](3,2)=(1 - gamma)*state[1]*state[2]/(state[0]*state[0]);
  (*m_A)[0](3,3)=gamma*state[1]/state[0];

//   (*m_A)[1][0][0]=0;
//   (*m_A)[1][0][1]=0;
  (*m_A)[1](0,2)=1;
//   (*m_A)[1][0][3]=0;

  (*m_A)[1](1,0)=-state[1]*state[2]/(state[0]*state[0]);
  (*m_A)[1](1,1)=state[2]/state[0];
  (*m_A)[1](1,2)=state[1]/state[0];
//   (*m_A)[1][1][3]=0;

  (*m_A)[1](2,0)=-(state[2]*state[2] - 1.0/2.0*temp)/(state[0]*state[0]);
  (*m_A)[1](2,1)=(1 - gamma)*state[1]/state[0];
  (*m_A)[1](2,2)=(3 - gamma)*state[2]/state[0];
  (*m_A)[1](2,3)=gamma - 1;

  (*m_A)[1](3,0)=state[2]/(state[0]*state[0])*(-gamma*state[3] + temp/state[0]);
  (*m_A)[1](3,1)=(1 - gamma)*state[1]*state[2]/(state[0]*state[0]);
  (*m_A)[1](3,2)=(gamma*state[3] - (temp/2.0 + (gamma - 1)*state[2]*state[2])/state[0])/state[0];
  (*m_A)[1](3,3)=gamma*state[2]/state[0];
  return;
}

//////////////////////////////////////////////////////////////////////////////

void StdBaseSolve::compute_Amatrix3D(Framework::State state, std::vector < RealMatrix > *m_A)
{
  CFreal gamma = 1.4;
  CFreal gamma1 = gamma - 1.0;
  CFreal r=state[0];
  cf_assert (r > 0);
  CFreal u=state[1]/r;
  CFreal v=state[2]/r;
  CFreal w=state[3]/r;
  CFreal v2 = u*u + v*v + w*w;

//  (*m_A)[0][0][0]=0;
  (*m_A)[0](0,1)=1.0;
//  (*m_A)[0][0][2]=0;
//  (*m_A)[0][0][3]=0;
//  (*m_A)[0][0][4]=0;

  (*m_A)[0](1,0)= gamma1*v2/2.0 - u*u;
  (*m_A)[0](1,1)=(3.0 - gamma)*u;
  (*m_A)[0](1,2)=-gamma1*v;
  (*m_A)[0](1,3)=-gamma1*w;
  (*m_A)[0](1,4)=gamma1;

  (*m_A)[0](2,0)=-u*v;
  (*m_A)[0](2,1)=v;
  (*m_A)[0](2,2)=u;
//  (*m_A)[0][2][3]=0;
//  (*m_A)[0][2][4]=0;

  (*m_A)[0](3,0)=-u*w;
  (*m_A)[0](3,1)=w;
//  (*m_A)[0][3][2]=0;
  (*m_A)[0](3,3)=u;
//  (*m_A)[0][3][4]=0;

  (*m_A)[0](4,0)=u*(gamma1*v2 - gamma*state[4]/r);
  (*m_A)[0](4,1)=gamma*state[4]/r - gamma1*(u*u + v2/2.0);
  (*m_A)[0](4,2)=-gamma1*u*v;
  (*m_A)[0](4,3)=-gamma1*u*w;
  (*m_A)[0](4,4)=gamma*u;

//  (*m_A)[1][0][0]=0;
//  (*m_A)[1][0][1]=0;
  (*m_A)[1](0,2)=1.0;
//  (*m_A)[1][0][3]=0;
//  (*m_A)[1][0][4]=0;

  (*m_A)[1](1,0)=-u*v;
  (*m_A)[1](1,1)=v;
  (*m_A)[1](1,2)=u;
//  (*m_A)[1][1][3]=0;
//  (*m_A)[1][1][4]=0;

  (*m_A)[1](2,0)= gamma1*v2/2.0 - v*v;
  (*m_A)[1](2,1)=-gamma1*u;
  (*m_A)[1](2,2)=(3.0 - gamma)*v;
  (*m_A)[1](2,3)=-gamma1*w;
  (*m_A)[1](2,4)=gamma1;

  (*m_A)[1](3,0)=-v*w;
//  (*m_A)[1][3][1]=0;
  (*m_A)[1](3,2)=w;
  (*m_A)[1](3,3)=v;
//  (*m_A)[1][3][4]=0;

  (*m_A)[1](4,0)=v*(gamma1*v2 - gamma*state[4]/r);
  (*m_A)[1](4,1)=-gamma1*u*v;
  (*m_A)[1](4,2)=gamma*state[4]/r - gamma1*(v*v + v2/2.0);
  (*m_A)[1](4,3)=-gamma1*v*w;
  (*m_A)[1](4,4)=gamma*v;

//  (*m_A)[2][0][0]=0;
//  (*m_A)[2][0][1]=0;
//  (*m_A)[2][0][2]=0;
  (*m_A)[2](0,3)=1.0;
//  (*m_A)[2][0][4]=0;

  (*m_A)[2](1,0)=-u*w;
  (*m_A)[2](1,1)=w;
//  (*m_A)[2][1][2]=0;
  (*m_A)[2](1,3)=u;
//  (*m_A)[2][1][4]=0;

  (*m_A)[2](2,0)=-v*w;
//  (*m_A)[2][2][1]=0;
  (*m_A)[2](2,2)=w;
  (*m_A)[2](2,3)=v;
//  (*m_A)[2][2][4]=0;

  (*m_A)[2](3,0)= gamma1*v2/2.0 - w*w;
  (*m_A)[2](3,1)=-gamma1*u;
  (*m_A)[2](3,2)=-gamma1*v;
  (*m_A)[2](3,3)=(3.0 - gamma)*w;
  (*m_A)[2](3,4)=gamma1;

  (*m_A)[2](4,0)=w*(gamma1*v2 - gamma*state[4]/r);
  (*m_A)[2](4,1)=-gamma1*u*w;
  (*m_A)[2](4,2)=-gamma1*v*w;
  (*m_A)[2](4,3)=gamma*state[4]/r - gamma1*(w*w + v2/2.0);
  (*m_A)[2](4,4)=gamma*w;

  return;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

