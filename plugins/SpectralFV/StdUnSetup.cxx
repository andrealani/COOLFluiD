#include "SpectralFV/StdUnSetup.hh"
#include "SpectralFV/SpectralFV.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFV {

Framework::MethodCommandProvider< StdUnSetup,SpectralFVMethodData,SpectralFVModule >
  stdUnSetupProvider("StdUnSetup");

//////////////////////////////////////////////////////////////////////////////

StdUnSetup::StdUnSetup(const std::string& name) :
  SpectralFVMethodCom(name),
  socket_faceNormTransfMatrices("faceNormTransfMatrices"),
  socket_volumes("volumes"),
  socket_faceSurf("faceSurf"),
  socket_gradients("gradients")
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
  StdUnSetup::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;

  result.push_back(&socket_faceNormTransfMatrices);
  result.push_back(&socket_volumes);
  result.push_back(&socket_faceSurf);
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

  // Force deallocate faceNormTransfMatrices
  DataHandle< RealMatrix > faceNormTransfMatrices = socket_faceNormTransfMatrices.getDataHandle();
  faceNormTransfMatrices.resize(0);

  // Force deallocate face surfaces
  DataHandle< CFreal > faceSurf = socket_faceSurf.getDataHandle();
  faceSurf.resize(0);

  // Force deallocate gradients
  DataHandle< vector< RealVector > > gradients = socket_gradients.getDataHandle();
  gradients.resize(0);

  // VARIABLES IN SPECTRALFVMETHODDATA

  // clear svLocalData
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();
  const CFuint nbrElemTypes = svLocalData.size();
  for (CFuint iElemType = 0; iElemType < nbrElemTypes; ++iElemType)
  {
    delete svLocalData[iElemType];
  }
  svLocalData.resize(0);

  // clear start index of inner faces with a certain orientation
  vector< CFuint >& innerFacesStartIdxs = getMethodData().getInnerFacesStartIdxs();
  innerFacesStartIdxs.resize(0);

  // clear start index of boundary faces with a certain orientation
  map<std::string, vector< vector< CFuint > > >& bndFacesStartIdxs = getMethodData().getBndFacesStartIdxs();
  bndFacesStartIdxs.clear();

}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV
}  // namespace COOLFluiD
