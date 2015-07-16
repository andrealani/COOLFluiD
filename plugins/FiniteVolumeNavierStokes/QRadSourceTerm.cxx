#include "FiniteVolumeNavierStokes/QRadSourceTerm.hh"
#include "NavierStokes/Euler2DVarSet.hh"
#include "Common/CFLog.hh"
#include "Framework/GeometricEntity.hh"
#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<QRadSourceTerm,
		       CellCenterFVMData, 
		       ComputeSourceTerm<CellCenterFVMData>, 
		       FiniteVolumeNavierStokesModule>
qradSourceFVMCCProvider("QRadST");
      
//////////////////////////////////////////////////////////////////////////////

void QRadSourceTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint > ("TID", "Temperature ID");
  options.addConfigOption< CFreal > ("RadRelaxationFactor", "Radiation Relaxation Factor");
}

//////////////////////////////////////////////////////////////////////////////

QRadSourceTerm::QRadSourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  socket_qrad("qrad")
{
  addConfigOptionsTo(this);
  
  m_TID = 0;
  setParameter("TID", &m_TID);

  m_relaxationFactor = 1.;
  setParameter("RadRelaxationFactor", &m_relaxationFactor);

}

//////////////////////////////////////////////////////////////////////////////

QRadSourceTerm::~QRadSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void QRadSourceTerm::setup()
{
  ComputeSourceTermFVMCC::setup();

  if (m_TID == 0) {
    m_TID = PhysicalModelStack::getActive()->getDim() + 1;
  }
}

//////////////////////////////////////////////////////////////////////////////

void QRadSourceTerm::computeSource(Framework::GeometricEntity *const element,
				   RealVector& source,
				   RealMatrix& jacobian)
{
  CFLogDebugMin( "QRadSourceTerm::computeSource()" << "\n");
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle<CFreal> qrad = socket_qrad.getDataHandle();
  const CFuint cellID = element->getID();
  State *const currState = element->getState(0);
  const CFreal r = (getMethodData().isAxisymmetric()) ? currState->getCoordinates()[YY] : 1.0;
  source[m_TID] += source[m_TID]*(1.-m_relaxationFactor) + m_relaxationFactor*qrad[cellID]*volumes[cellID]*r;
  CFLog(DEBUG_MAX, "QRadSourceTerm::computeSource() => source = " << source << "\n");
  CFLog(DEBUG_MAX, "QRadSourceTerm::computeSource() => m_relaxationFactor  = " << m_relaxationFactor << "\n");

}

//////////////////////////////////////////////////////////////////////////////
 
vector<SafePtr<BaseDataSocketSink> > QRadSourceTerm::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = ComputeSourceTermFVMCC::needsSockets();
  result.push_back(&socket_qrad);
  return result;
}

//////////////////////////////////////////////////////////////////////////////
  
    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
