#include "ConstantPolyRec.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<ConstantPolyRec, 
		       CellCenterFVMData,
		       PolyReconstructor<CellCenterFVMData>,
		       FiniteVolumeModule> 
constantPolyRecProvider("Constant");

//////////////////////////////////////////////////////////////////////////////

ConstantPolyRec::ConstantPolyRec(const std::string& name) :
  FVMCC_PolyRec(name)
{
}

//////////////////////////////////////////////////////////////////////////////

ConstantPolyRec::~ConstantPolyRec()
{
}
      
//////////////////////////////////////////////////////////////////////////////

void ConstantPolyRec::extrapolateImpl(GeometricEntity* const face)
{
  FVMCC_PolyRec::baseExtrapolateImpl(face);
  FVMCC_PolyRec::constantRecImpl(face);
}

//////////////////////////////////////////////////////////////////////////////

void ConstantPolyRec::extrapolateImpl(GeometricEntity* const face,
				      CFuint iVar, CFuint leftOrRight)
{
  FVMCC_PolyRec::constantRecImpl(face, iVar, leftOrRight);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
