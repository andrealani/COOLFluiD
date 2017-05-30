#include "FiniteVolume/RMSJouleHeatSourceTerm.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/CellCenterFVMData.hh"
#include "FiniteVolume/FiniteVolume.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<RMSJouleHeatSourceTerm,
		       CellCenterFVMData, 
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeModule>
rmsJouleHeatSTFVMCCProvider("RMSJouleHeatST");

//////////////////////////////////////////////////////////////////////////////

RMSJouleHeatSourceTerm::
RMSJouleHeatSourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  socket_rmsJouleHeatSource("rmsJouleHeatSource")
{
}

//////////////////////////////////////////////////////////////////////////////

RMSJouleHeatSourceTerm::~RMSJouleHeatSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
RMSJouleHeatSourceTerm::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = ComputeSourceTermFVMCC::needsSockets();

  result.push_back(&socket_rmsJouleHeatSource);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void RMSJouleHeatSourceTerm::setup()
{
  using namespace COOLFluiD::Framework;

  ComputeSourceTermFVMCC::setup();

}

//////////////////////////////////////////////////////////////////////////////

void RMSJouleHeatSourceTerm::computeSource
(Framework::GeometricEntity *const element, RealVector& source, RealMatrix& jacobian)
{
  using namespace COOLFluiD::Framework;

  const CFuint nbEqs =
    PhysicalModelStack::getActive()->getEquationSubSysDescriptor().getNbEqsSS();

  // this is needed for coupling
  if (nbEqs == 4 || nbEqs == source.size()) {
    CFLogDebugMin( "RMSJouleHeatSourceTerm::computeSource()" << "\n");

    const CFuint elemID = element->getID();

    DataHandle<CFreal> rmsJouleHeatSource = socket_rmsJouleHeatSource.getDataHandle();
    DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

    // is not perturbed because it is computed in command, here is got just data handle
    source[3] = rmsJouleHeatSource[elemID]*volumes[elemID];

    #ifdef DEBUG
      CFout <<"source * volumes (in RMSJouleHeatSourceTerm) = " << source << "\n"; 
    #endif
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
