#include "Common/PE.hh"
#include "Common/BadValueException.hh"

#include "MathTools/MathConsts.hh"

#include "Environment/FileHandlerOutput.hh"
#include "Environment/DirPaths.hh"
#include "Environment/SingleBehaviorFactory.hh"

#include "Framework/PathAppender.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalConsts.hh"

#include "FiniteVolume/CellCenterFVM.hh"
#include "FiniteVolume/SphericalDerivatives.hh"
#include "FiniteVolume/FiniteVolume.hh"

/////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Numerics::FiniteVolume;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////
      
MethodCommandProvider<SphericalDerivatives, DataProcessingData, FiniteVolumeModule>
sphericalDerivativesProvider("SphericalDerivatives");

//////////////////////////////////////////////////////////////////////////////

void SphericalDerivatives::defineConfigOptions(Config::OptionList& options)
{
  //options.addConfigOption< CFuint >("SaveRate","Save Output File every...iterations.");
}

//////////////////////////////////////////////////////////////////////////////

SphericalDerivatives::SphericalDerivatives(const std::string& name) :
  DataProcessingCom(name),
  socket_uR("uR"),
  socket_uTheta("uTheta"),
  socket_uPhi("uPhi"),
  socket_uX("uX"),
  socket_uY("uY"),
  socket_uZ("uZ"),
  socket_states("states")
{
  addConfigOptionsTo(this);
  
  // m_saveRate = 1;
  // setParameter("SaveRate",&m_saveRate);
}

//////////////////////////////////////////////////////////////////////////////

SphericalDerivatives::~SphericalDerivatives()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
SphericalDerivatives::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  result.push_back(&socket_uR);
  result.push_back(&socket_uTheta);
  result.push_back(&socket_uPhi);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
SphericalDerivatives::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_states);
  result.push_back(&socket_uX);
  result.push_back(&socket_uY);
  result.push_back(&socket_uZ);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void SphericalDerivatives::setup()
{
  CFAUTOTRACE;
  
  const CFuint usize = socket_uX.getDataHandle().size();
  socket_uR.getDataHandle().resize(usize);
  socket_uTheta.getDataHandle().resize(usize);
  socket_uPhi.getDataHandle().resize(usize);
}

//////////////////////////////////////////////////////////////////////////////

void SphericalDerivatives::execute()
{
  CFLog(VERBOSE, "SphericalDerivatives::execute() => START\n");
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<CFreal> uR = socket_uR.getDataHandle();
  DataHandle<CFreal> uTheta = socket_uTheta.getDataHandle();
  DataHandle<CFreal> uPhi = socket_uPhi.getDataHandle();
  DataHandle<CFreal> uX = socket_uX.getDataHandle();
  DataHandle<CFreal> uY = socket_uY.getDataHandle();
  DataHandle<CFreal> uZ = socket_uZ.getDataHandle();

  const CFuint nbCells = states.size();
  const CFuint nbVars = uX.size()/nbCells;
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  for (CFuint i = 0; i < nbCells; ++i) {
    const RealVector& coord = states[i]->getCoordinates();
    const CFreal x = coord[XX];
    const CFreal y = coord[YY];
    const CFreal r = std::sqrt(x*x + y*y);
    if (dim == DIM_2D) {
      const CFuint start = i*nbVars;
      for (CFuint iVar = 0; iVar < nbVars; ++iVar) {
	const CFuint idx = start + iVar;
	uR[idx]     = x/r*uX[idx] + y/r*uY[idx];
	uTheta[idx] = -y*uX[idx] + x*uY[idx];
	uPhi[idx]   = 0.;
      }
    }
    
    if (dim == DIM_3D) {
      const CFreal z = coord[ZZ];
      const CFreal r3 = std::sqrt(x*x + y*y + z*z);
      const CFuint start = i*nbVars;
      for (CFuint iVar = 0; iVar < nbVars; ++iVar) {
	const CFuint idx = start + iVar;
	uR[idx]     = x/r3*uX[idx] + y/r3*uY[idx] + z/r3*uZ[idx];
	uTheta[idx] = -y*uX[idx] + x*uY[idx];
	uPhi[idx]   = z*x/r*uX[idx] + z*y/r*uY[idx] - r*uZ[idx];
      }
    }
  }

  CFLog(VERBOSE, "SphericalDerivatives::execute() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

void SphericalDerivatives::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

