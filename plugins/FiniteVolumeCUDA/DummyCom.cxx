#include "Framework/MethodCommandProvider.hh"
#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolumeCUDA/DummyCom.hh"
#include "FiniteVolumeCUDA/FiniteVolumeCUDAParalution.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<DummyCom, CellCenterFVMData,FiniteVolumeCUDAParalutionModule> 
dummyComProvider("DummyCom");

//////////////////////////////////////////////////////////////////////////////

DummyCom::DummyCom(const std::string& name) :
  CellCenterFVMCom(name)
{
}
 
//////////////////////////////////////////////////////////////////////////////

DummyCom::~DummyCom()
{
}

//////////////////////////////////////////////////////////////////////////////

void DummyCom::execute()
{
  CFAUTOTRACE;
  
  CFLog(VERBOSE, "DummyCom::execute() => start\n");
  CFLog(VERBOSE, "DummyCom::execute() => end\n");
}
      
//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
