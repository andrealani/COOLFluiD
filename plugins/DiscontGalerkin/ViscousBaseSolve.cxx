#include "Common/COOLFluiD.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"

#include "DiscontGalerkin/ViscousBaseSolve.hh"
#include "DiscontGalerkin/DiscontGalerkin.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ViscousBaseSolve,DiscontGalerkinSolverData,DiscontGalerkinModule >
  viscousBaseSolveProvider("ViscousBaseSolve");

//////////////////////////////////////////////////////////////////////////////

ViscousBaseSolve::ViscousBaseSolve(const std::string& name)
  : StdBaseSolve(name)
{
}

//////////////////////////////////////////////////////////////////////////////

ViscousBaseSolve::~ViscousBaseSolve()
{
}

//////////////////////////////////////////////////////////////////////////////

void ViscousBaseSolve::setup()
{
  CFAUTOTRACE;
  StdBaseSolve::setup();
  m_Re = getMethodData().getReynolds();
}

//////////////////////////////////////////////////////////////////////////////

void ViscousBaseSolve::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void ViscousBaseSolve::compute_Kmatrix2D(Framework::State state, std::vector < std::vector < RealMatrix > > *m_K)
{
  CFreal Pr = 0.72;
  CFreal kappa = 1.4;

  CFreal ro=state[0];
  cf_assert(ro > 0.);
  CFreal v1=state[1]/state[0];
  CFreal v2=state[2]/state[0];

  CFreal oneOverm_Rew1 = 1.0/(m_Re*ro);

  CFreal hlpK51 = -(v1*v1+v2*v2)*oneOverm_Rew1
                  + kappa*oneOverm_Rew1/Pr*(-state[3]/ro + v1*v1 + v2*v2);
  CFreal hlpK52 = (1.-kappa/Pr)*v1*oneOverm_Rew1;
  CFreal hlpK53 = (1.-kappa/Pr)*v2*oneOverm_Rew1;
  CFreal hlpK55 =  kappa/Pr*oneOverm_Rew1;

  (*m_K)[0][0](1,0)=-4.0/3.0*v1*oneOverm_Rew1;
  (*m_K)[0][0](1,1)=4.0/3.0*oneOverm_Rew1;
  (*m_K)[0][0](2,0)=-v2*oneOverm_Rew1;
  (*m_K)[0][0](2,2)=oneOverm_Rew1;
  (*m_K)[0][0](3,0)=hlpK51-1.0/3.0*v1*v1*oneOverm_Rew1;
  (*m_K)[0][0](3,1)=hlpK52 + 1.0/3.0*v1*oneOverm_Rew1;
  (*m_K)[0][0](3,2)=hlpK53;
  (*m_K)[0][0](3,3)=hlpK55;

  (*m_K)[1][1](1,0)=-v1*oneOverm_Rew1;
  (*m_K)[1][1](1,1)= oneOverm_Rew1;
  (*m_K)[1][1](2,0)=-4.0/3.0*v2*oneOverm_Rew1;
  (*m_K)[1][1](2,2)=4.0/3.0*oneOverm_Rew1;
  (*m_K)[1][1](3,0)=hlpK51-1.0/3.0*v2*v2*oneOverm_Rew1;
  (*m_K)[1][1](3,1)=hlpK52;
  (*m_K)[1][1](3,2)=hlpK53 + 1.0/3.0*v2*oneOverm_Rew1;
  (*m_K)[1][1](3,3)=hlpK55;

  (*m_K)[0][1](1,0)= 2.0/3.0*v2*oneOverm_Rew1;
  (*m_K)[0][1](1,2)=-2.0/3.0*oneOverm_Rew1;
  (*m_K)[0][1](2,0)=-v1*oneOverm_Rew1;
  (*m_K)[0][1](2,1)= oneOverm_Rew1;
  (*m_K)[0][1](3,0)=-1.0/3.0*v1*v2*oneOverm_Rew1;
  (*m_K)[0][1](3,1)= v2*oneOverm_Rew1;
  (*m_K)[0][1](3,2)=-2.0/3.0*v1*oneOverm_Rew1;

  (*m_K)[1][0](1,0)=-v2*oneOverm_Rew1;
  (*m_K)[1][0](1,2)= oneOverm_Rew1;
  (*m_K)[1][0](2,0)= 2.0/3.0*v1*oneOverm_Rew1;
  (*m_K)[1][0](2,1)=-2.0/3.0*oneOverm_Rew1;
  (*m_K)[1][0](3,0)=-1.0/3.0*v1*v2*oneOverm_Rew1;
  (*m_K)[1][0](3,1)=-2.0/3.0*v2*oneOverm_Rew1;
  (*m_K)[1][0](3,2)= v1*oneOverm_Rew1;

  return;
}

//////////////////////////////////////////////////////////////////////////////

void ViscousBaseSolve::compute_Kmatrix3D(Framework::State state,  std::vector < std::vector < RealMatrix > > *m_K)
{
  CFreal Pr = 0.72;
  CFreal kappa = 1.4;

  CFreal ro=state[0];
  cf_assert(ro > 0.);
  CFreal v1=state[1]/state[0];
  CFreal v2=state[2]/state[0];
  CFreal v3=state[3]/state[0];

  CFreal oneOverm_Rew1 = 1.0/(m_Re*ro);

  CFreal hlpK51 = -(v1*v1+v2*v2+v3*v3)*oneOverm_Rew1
                  + kappa*oneOverm_Rew1/Pr*(-state[4]/ro + v1*v1 + v2*v2 + v3*v3);
  CFreal hlpK52 = (1.-kappa/Pr)*v1*oneOverm_Rew1;
  CFreal hlpK53 = (1.-kappa/Pr)*v2*oneOverm_Rew1;
  CFreal hlpK54 = (1.-kappa/Pr)*v3*oneOverm_Rew1;
  CFreal hlpK55 =  kappa/Pr*oneOverm_Rew1;

  (*m_K)[0][0](1,0)=-4.0/3.0*v1*oneOverm_Rew1;
  (*m_K)[0][0](1,1)=4.0/3.0*oneOverm_Rew1;
  (*m_K)[0][0](2,0)=-v2*oneOverm_Rew1;
  (*m_K)[0][0](2,2)=oneOverm_Rew1;
  (*m_K)[0][0](3,0)=-v3*oneOverm_Rew1;
  (*m_K)[0][0](3,3)= oneOverm_Rew1;
  (*m_K)[0][0](4,0)=hlpK51-1.0/3.0*v1*v1*oneOverm_Rew1;
  (*m_K)[0][0](4,1)=hlpK52 + 1.0/3.0*v1*oneOverm_Rew1;
  (*m_K)[0][0](4,2)=hlpK53;
  (*m_K)[0][0](4,3)=hlpK54;
  (*m_K)[0][0](4,4)=hlpK55;

  (*m_K)[1][1](1,0)=-v1*oneOverm_Rew1;
  (*m_K)[1][1](1,1)= oneOverm_Rew1;
  (*m_K)[1][1](2,0)=-4.0/3.0*v2*oneOverm_Rew1;
  (*m_K)[1][1](2,2)=4.0/3.0*oneOverm_Rew1;
  (*m_K)[1][1](3,0)=-v3*oneOverm_Rew1;
  (*m_K)[1][1](3,3)=oneOverm_Rew1;
  (*m_K)[1][1](4,0)=hlpK51-1.0/3.0*v2*v2*oneOverm_Rew1;
  (*m_K)[1][1](4,1)=hlpK52;
  (*m_K)[1][1](4,2)=hlpK53 + 1.0/3.0*v2*oneOverm_Rew1;
  (*m_K)[1][1](4,3)=hlpK54;
  (*m_K)[1][1](4,4)=hlpK55;

  (*m_K)[2][2](1,0)=-v1*oneOverm_Rew1;
  (*m_K)[2][2](1,1)= oneOverm_Rew1;
  (*m_K)[2][2](2,0)=-v2*oneOverm_Rew1;
  (*m_K)[2][2](2,2)= oneOverm_Rew1;
  (*m_K)[2][2](3,0)=-4.0/3.0*v3*oneOverm_Rew1;
  (*m_K)[2][2](3,3)= 4.0/3.0*oneOverm_Rew1;
  (*m_K)[2][2](4,0)=hlpK51-1.0/3.0*v3*v3*oneOverm_Rew1;
  (*m_K)[2][2](4,1)=hlpK52;
  (*m_K)[2][2](4,2)=hlpK53;
  (*m_K)[2][2](4,3)=hlpK54 + 1.0/3.0*v3*oneOverm_Rew1;
  (*m_K)[2][2](4,4)=hlpK55;

  (*m_K)[0][1](1,0)= 2.0/3.0*v2*oneOverm_Rew1;
  (*m_K)[0][1](1,2)=-2.0/3.0*oneOverm_Rew1;
  (*m_K)[0][1](2,0)=-v1*oneOverm_Rew1;
  (*m_K)[0][1](2,1)= oneOverm_Rew1;
  (*m_K)[0][1](4,0)=-1.0/3.0*v1*v2*oneOverm_Rew1;
  (*m_K)[0][1](4,1)= v2*oneOverm_Rew1;
  (*m_K)[0][1](4,2)=-2.0/3.0*v1*oneOverm_Rew1;

  (*m_K)[0][2](1,0)= 2.0/3.0*v3*oneOverm_Rew1;
  (*m_K)[0][2](1,2)=-2.0/3.0*oneOverm_Rew1;
  (*m_K)[0][2](3,0)=-v1*oneOverm_Rew1;
  (*m_K)[0][2](3,1)= oneOverm_Rew1;
  (*m_K)[0][2](4,0)=-1.0/3.0*v1*v3*oneOverm_Rew1;
  (*m_K)[0][2](4,1)= v3*oneOverm_Rew1;
  (*m_K)[0][2](4,3)=-2.0/3.0*v1*oneOverm_Rew1;

  (*m_K)[1][0](1,0)=-v2*oneOverm_Rew1;
  (*m_K)[1][0](1,2)= oneOverm_Rew1;
  (*m_K)[1][0](2,0)= 2.0/3.0*v1*oneOverm_Rew1;
  (*m_K)[1][0](2,1)=-2.0/3.0*oneOverm_Rew1;
  (*m_K)[1][0](4,0)=-1.0/3.0*v1*v2*oneOverm_Rew1;
  (*m_K)[1][0](4,1)=-2.0/3.0*v2*oneOverm_Rew1;
  (*m_K)[1][0](4,2)= v1*oneOverm_Rew1;

  (*m_K)[1][2](2,0)= 2.0/3.0*v3*oneOverm_Rew1;
  (*m_K)[1][2](2,2)=-2.0/3.0*oneOverm_Rew1;
  (*m_K)[1][2](3,0)=-v2*oneOverm_Rew1;
  (*m_K)[1][2](3,2)= oneOverm_Rew1;
  (*m_K)[1][2](4,0)=-1.0/3.0*v2*v3*oneOverm_Rew1;
  (*m_K)[1][2](4,2)= v3*oneOverm_Rew1;
  (*m_K)[1][2](4,3)=-2.0/3.0*v2*oneOverm_Rew1;

  (*m_K)[2][0](1,0)=-v3*oneOverm_Rew1;
  (*m_K)[2][0](1,3)= oneOverm_Rew1;
  (*m_K)[2][0](3,0)= 2.0/3.0*v1*oneOverm_Rew1;
  (*m_K)[2][0](3,1)=-2.0/3.0*oneOverm_Rew1;
  (*m_K)[2][0](4,0)=-1.0/3.0*v1*v3*oneOverm_Rew1;
  (*m_K)[2][0](4,1)=-2.0/3.0*v3*oneOverm_Rew1;
  (*m_K)[2][0](4,3)= v1*oneOverm_Rew1;

  (*m_K)[2][1](2,0)=-v3*oneOverm_Rew1;
  (*m_K)[2][1](2,3)= oneOverm_Rew1;
  (*m_K)[2][1](3,0)= 2.0/3.0*v2*oneOverm_Rew1;
  (*m_K)[2][1](3,2)=-2.0/3.0*oneOverm_Rew1;
  (*m_K)[2][1](4,0)=-1.0/3.0*v2*v3*oneOverm_Rew1;
  (*m_K)[2][1](4,2)=-2.0/3.0*v3*oneOverm_Rew1;
  (*m_K)[2][1](4,3)= v2*oneOverm_Rew1;

  return;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

