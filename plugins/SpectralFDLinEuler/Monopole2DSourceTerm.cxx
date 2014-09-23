#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

#include "SpectralFDLinEuler/Monopole2DSourceTerm.hh"
#include "SpectralFDLinEuler/SpectralFDLinEuler.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<Monopole2DSourceTerm, SpectralFDMethodData, SpectralFDLinEulerModule>
 Monopole2DSourceTermProvider("Monopole2DSourceTerm");

//////////////////////////////////////////////////////////////////////////////

Monopole2DSourceTerm::Monopole2DSourceTerm(const std::string& name) :
  StdSourceTerm(name),
  m_currTime()
{
}

//////////////////////////////////////////////////////////////////////////////

Monopole2DSourceTerm::~Monopole2DSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void Monopole2DSourceTerm::getSourceTermData()
{
  StdSourceTerm::getSourceTermData();

  // get current time
  m_currTime = SubSystemStatusStack::getActive()->getCurrentTimeDim();
}

//////////////////////////////////////////////////////////////////////////////

void Monopole2DSourceTerm::addSourceTerm()
{
  // pre.minary data
  const CFreal alpha = log(2.0)/2.0;
  const CFreal eps1 = 0.5*sqrt(1.4);
  const CFreal eps2 = 0.5*1.4*sqrt(1.4);
  const CFreal om = MathTools::MathConsts::CFrealPi()/15.0*sqrt(1.4);
  const CFreal timeFac1 = eps1*sin(om*m_currTime);
  const CFreal timeFac2 = eps2*sin(om*m_currTime);

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
    rhs[resID  ] += resFactor*m_solPntJacobDets[iSol]*exp(-alpha*(x*x+y*y))*timeFac1;
    rhs[resID+3] += resFactor*m_solPntJacobDets[iSol]*exp(-alpha*(x*x+y*y))*timeFac2;
  }
}

//////////////////////////////////////////////////////////////////////////////

void Monopole2DSourceTerm::setup()
{
  CFAUTOTRACE;
  StdSourceTerm::setup();
}

//////////////////////////////////////////////////////////////////////////////

void Monopole2DSourceTerm::unsetup()
{
  CFAUTOTRACE;
  StdSourceTerm::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    Monopole2DSourceTerm::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = StdSourceTerm::needsSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
