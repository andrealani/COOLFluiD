#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

#include "SpectralFDLinEuler/Dipole2DSourceTerm.hh"
#include "SpectralFDLinEuler/SpectralFDLinEuler.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<Dipole2DSourceTerm, SpectralFDMethodData, SpectralFDLinEulerModule>
 Dipole2DSourceTermProvider("Dipole2DSourceTerm");

//////////////////////////////////////////////////////////////////////////////

Dipole2DSourceTerm::Dipole2DSourceTerm(const std::string& name) :
  StdSourceTerm(name),
  m_currTime()
{
}

//////////////////////////////////////////////////////////////////////////////

Dipole2DSourceTerm::~Dipole2DSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void Dipole2DSourceTerm::getSourceTermData()
{
  StdSourceTerm::getSourceTermData();

  // get current time
  m_currTime = SubSystemStatusStack::getActive()->getCurrentTimeDim();
}

//////////////////////////////////////////////////////////////////////////////

void Dipole2DSourceTerm::addSourceTerm()
{
  // pre.minary data
  const CFreal alpha = log(2.0)/5.0;
  const CFreal beta = MathTools::MathConsts::CFrealPi()/10.0;
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
    if (std::abs(x) <= 5.0)
    {
      const CFreal y = ((*m_cellStates)[iSol]->getCoordinates())[YY]*refLength;
      rhs[resID+1] += resFactor*m_solPntJacobDets[iSol]*exp(-alpha*y*y)*cos(beta*x)*timeFac;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Dipole2DSourceTerm::setup()
{
  CFAUTOTRACE;
  StdSourceTerm::setup();
}

//////////////////////////////////////////////////////////////////////////////

void Dipole2DSourceTerm::unsetup()
{
  CFAUTOTRACE;
  StdSourceTerm::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    Dipole2DSourceTerm::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = StdSourceTerm::needsSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
