#include "Framework/MethodCommandProvider.hh"

#include "SpectralFD/StdUnSetup.hh"
#include "SpectralFD/SpectralFD.hh"
#include "SpectralFD/SpectralFDElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFD {

Framework::MethodCommandProvider< StdUnSetup,SpectralFDMethodData,SpectralFDModule >
  stdUnSetupProvider("StdUnSetup");

//////////////////////////////////////////////////////////////////////////////

StdUnSetup::StdUnSetup(const std::string& name) :
  SpectralFDMethodCom(name),
  socket_volumes("volumes"),
  socket_faceJacobVecSizeFaceFlxPnts("faceJacobVecSizeFaceFlxPnts"),
  socket_gradients("gradients")
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
  StdUnSetup::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;

  result.push_back(&socket_volumes);
  result.push_back(&socket_faceJacobVecSizeFaceFlxPnts);
  result.push_back(&socket_gradients);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  StdUnSetup::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result;
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdUnSetup::execute()
{
  CFAUTOTRACE;

  // SOCKETS

  // Force deallocate volumes
  DataHandle< CFreal > volumes = socket_volumes.getDataHandle();
  volumes.resize(0);

  // Force deallocate socket_faceJacobVecSizeFaceFlxPnts
  DataHandle< vector< CFreal > > faceJacobVecSizeFaceFlxPnts = socket_faceJacobVecSizeFaceFlxPnts.getDataHandle();
  faceJacobVecSizeFaceFlxPnts.resize(0);

  // Force deallocate gradients
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();
  gradients.resize(0);

  // VARIABLES IN SPECTRALFDMETHODDATA

  // clear sdLocalData
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();
  const CFuint nbrElemTypes = sdLocalData.size();
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    delete sdLocalData[iElemType];
  }
  sdLocalData.resize(0);

  // clear start index of inner faces with a certain orientation
  vector< CFuint >& innerFacesStartIdxs = getMethodData().getInnerFacesStartIdxs();
  innerFacesStartIdxs.resize(0);

  // clear start index of boundary faces with a certain orientation
  map<std::string, vector< vector< CFuint > > >& bndFacesStartIdxs = getMethodData().getBndFacesStartIdxs();
  bndFacesStartIdxs.clear();

}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD
}  // namespace COOLFluiD
