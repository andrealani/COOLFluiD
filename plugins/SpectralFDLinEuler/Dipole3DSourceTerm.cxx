#include "Common/CFLog.hh"

#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SubSystemStatus.hh"

#include "SpectralFDLinEuler/Dipole3DSourceTerm.hh"
#include "SpectralFDLinEuler/SpectralFDLinEuler.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<Dipole3DSourceTerm, SpectralFDMethodData, SpectralFDLinEulerModule>
 Dipole3DSourceTermProvider("Dipole3DSourceTerm");

//////////////////////////////////////////////////////////////////////////////

Dipole3DSourceTerm::Dipole3DSourceTerm(const std::string& name) :
  StdSourceTerm(name),
  m_currTime()
{
}

//////////////////////////////////////////////////////////////////////////////

Dipole3DSourceTerm::~Dipole3DSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void Dipole3DSourceTerm::getSourceTermData()
{
  StdSourceTerm::getSourceTermData();

  // get current time
  m_currTime = SubSystemStatusStack::getActive()->getCurrentTimeDim();
}

//////////////////////////////////////////////////////////////////////////////

void Dipole3DSourceTerm::addSourceTerm()
{
  // pre.minary data
  const CFreal alpha = log(2.)/0.0125;
  const CFreal eps = 0.001;
  const CFreal om = MathTools::MathConsts::CFrealPi();

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
    const CFreal x = ((*m_cellStates)[iSol]->getCoordinates())[XX]*refLength-0.5;
    const CFreal y = ((*m_cellStates)[iSol]->getCoordinates())[YY]*refLength-0.5;
    const CFreal z = ((*m_cellStates)[iSol]->getCoordinates())[ZZ]*refLength-0.5;
    const CFreal f = exp(-alpha*(x*x+y*y+z*z));

    rhs[resID+1] += resFactor*m_solPntJacobDets[iSol]*(eps*f*cos(MathTools::MathConsts::CFrealPi()*x)*sin(om*m_currTime));
  }
}

//////////////////////////////////////////////////////////////////////////////

void Dipole3DSourceTerm::setup()
{
  CFAUTOTRACE;
  StdSourceTerm::setup();
}

//////////////////////////////////////////////////////////////////////////////

void Dipole3DSourceTerm::unsetup()
{
  CFAUTOTRACE;
  StdSourceTerm::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
    Dipole3DSourceTerm::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = StdSourceTerm::needsSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD
