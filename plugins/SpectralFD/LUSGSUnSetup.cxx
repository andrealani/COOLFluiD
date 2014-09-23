#include "SpectralFD/LUSGSUnSetup.hh"
#include "SpectralFD/SpectralFD.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFD {

Framework::MethodCommandProvider< LUSGSUnSetup,SpectralFDMethodData,SpectralFDModule >
  lusgsUnSetupProvider("LUSGSUnSetup");

//////////////////////////////////////////////////////////////////////////////

LUSGSUnSetup::LUSGSUnSetup(const std::string& name) :
  StdUnSetup(name),
  socket_diagBlockJacobMatr("diagBlockJacobMatr"),
  socket_rhsCurrStatesSet("rhsCurrStatesSet"),
  socket_statesSetStateIDs("statesSetStateIDs")
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
  LUSGSUnSetup::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result = StdUnSetup::needsSockets();

  result.push_back(&socket_diagBlockJacobMatr     );
  result.push_back(&socket_rhsCurrStatesSet       );
  result.push_back(&socket_statesSetStateIDs      );

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
  LUSGSUnSetup::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result = StdUnSetup::providesSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void LUSGSUnSetup::execute()
{
  CFAUTOTRACE;

  // call execute() of the parent class
  StdUnSetup::execute();

  // SOCKETS

  // Force deallocate diagBlockJacobMatr
  DataHandle< RealMatrix > diagBlockJacobMatr = socket_diagBlockJacobMatr.getDataHandle();
  diagBlockJacobMatr.resize(0);

  // Force deallocate socket_rhsCurrStatesSet
  DataHandle< CFreal > rhsCurrStatesSet = socket_rhsCurrStatesSet.getDataHandle();
  rhsCurrStatesSet.resize(0);

  // Force deallocate statesSetStateIDs
  DataHandle< vector< CFuint > > statesSetStateIDs = socket_statesSetStateIDs.getDataHandle();
  statesSetStateIDs.resize(0);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD
}  // namespace COOLFluiD
