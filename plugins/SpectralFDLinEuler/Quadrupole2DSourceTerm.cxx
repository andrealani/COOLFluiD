#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

#include "SpectralFDLinEuler/Quadrupole2DSourceTerm.hh"
#include "SpectralFDLinEuler/SpectralFDLinEuler.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<Quadrupole2DSourceTerm, SpectralFDMethodData, SpectralFDLinEulerModule>
 Quadrupole2DSourceTermProvider("Quadrupole2DSourceTerm");

//////////////////////////////////////////////////////////////////////////////

Quadrupole2DSourceTerm::Quadrupole2DSourceTerm(const std::string& name) :
  StdSourceTerm(name),
  m_currTime()
{
}

//////////////////////////////////////////////////////////////////////////////

Quadrupole2DSourceTerm::~Quadrupole2DSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void Quadrupole2DSourceTerm::getSourceTermData()
{
  StdSourceTerm::getSourceTermData();

  // get current time
  m_currTime = SubSystemStatusStack::getActive()->getCurrentTimeDim();
}

//////////////////////////////////////////////////////////////////////////////

void Quadrupole2DSourceTerm::addSourceTerm()
{
  // pre.minary data
  const CFreal alpha = log(2)/5.0;
  const CFreal beta = MathTools::MathConsts::CFrealPi()/20.0;
  const CFreal eps = 0.014;
  const CFreal om = MathTools::MathConsts::CFrealPi()/30.0*sqrt(1.4);
  const CFreal timeFac = eps*sin(om*m_currTime);

  // get the datahandle of the rhs
  DataHandle< CFreal > rhs = socket_rhs.getDataHandle();

  // get residual factor
  const CFreal resFactor = getMethodData().getResFactor();

  // get reference length
  const CFreal refLength = PhysicalModelStack::getActive()->getImplementor()->getRefLength();

  // loop over solution points in this cell to add the dipole source term
  CFuint resID = m_nbrEqs*( (*m_cellStates)[0]->getLocalID() );
  const CFuint nbrSol = m_cellStates->size();
  for (CFuint iSol = 0; iSol < nbrSol; ++iSol, resID += m_nbrEqs)
  {
    const CFreal x = ((*m_cellStates)[iSol]->getCoordinates())[XX]*refLength;
    const CFreal y = ((*m_cellStates)[iSol]->getCoordinates())[YY]*refLength;
    if (std::abs(x) <= 10.0 && std::abs(y) <= 10.0)
    {
      rhs[resID+1] += resFactor*m_solPntJacobDets[iSol]*sin(beta*x)*exp(-alpha*y*y)*timeFac;
      rhs[resID+2] -= resFactor*m_solPntJacobDets[iSol]*sin(beta*y)*exp(-alpha*x*x)*timeFac;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Quadrupole2DSourceTerm::setup()
{
  CFAUTOTRACE;
  StdSourceTerm::setup();
}

//////////////////////////////////////////////////////////////////////////////

void Quadrupole2DSourceTerm::unsetup()
{
  CFAUTOTRACE;
  StdSourceTerm::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr< BaseDataSocketSink > >
    Quadrupole2DSourceTerm::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = StdSourceTerm::needsSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
