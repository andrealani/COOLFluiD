#include "ConstantSourceTerm.hh"
#include "Common/CFLog.hh"
#include "Framework/GeometricEntity.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////

MethodStrategyProvider<ConstantSourceTerm,
		      CellCenterFVMData, 
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeModule>
constantSTFVMCCProvider("ConstantST");

//////////////////////////////////////////////////////////////////////

ConstantSourceTerm::ConstantSourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name)
{
}

//////////////////////////////////////////////////////////////////////

ConstantSourceTerm::~ConstantSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////

void ConstantSourceTerm::setup()
{
  ComputeSourceTermFVMCC::setup();
}

//////////////////////////////////////////////////////////////////////

void ConstantSourceTerm::computeSource(GeometricEntity *const element,
				       RealVector& source,
				       RealMatrix& jacobian)
{
  CFLogDebugMin( "ConstantSourceTerm::computeSource()" << "\n");
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

  ///@todo make the constant an option
  source[0] = 1.*volumes[element->getID()];
}

//////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////
