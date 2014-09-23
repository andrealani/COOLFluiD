#include <fstream>
#include <iomanip>
#include <stdio.h>

#include "Common/COOLFluiD.hh"
#if CF_HAVE_UNISTD_H
  extern "C"
  {
    #include "TecplotWriter/TECXXX.h"
  }
#endif

#include "Common/PE.hh"
#include "Common/CFMap.hh"
#include "Common/BadValueException.hh"
#include "Common/OSystem.hh"

#include "Environment/DirPaths.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"

#include "Framework/ConvectiveVarSet.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/MapGeoEnt.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/StdTrsGeoBuilder.hh"
#include "Framework/DataHandleOutput.hh"
#include "Framework/PathAppender.hh"
#include "Framework/SubSystemStatus.hh"

#include "NavierStokes/EulerTerm.hh"
#include "NavierStokes/EulerVarSet.hh"

#include "TecplotWriterNavierStokes/TecplotWriterNavierStokes.hh"
#include "TecplotWriterNavierStokes/WriteInstantAndAvgSolutionHighOrder.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace boost::filesystem;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WriteInstantAndAvgSolutionHighOrder, TecWriterData, TecplotWriterNavierStokesModule>
writeInstantAndAvgSolutionHighOrderProvider("WriteInstantAndAvgSolutionHighOrder");

//////////////////////////////////////////////////////////////////////////////

void WriteInstantAndAvgSolutionHighOrder::defineConfigOptions(Config::OptionList& options)
{

}

//////////////////////////////////////////////////////////////////////////////

WriteInstantAndAvgSolutionHighOrder::WriteInstantAndAvgSolutionHighOrder(const std::string& name) : WriteSolution(name)
{

}

//////////////////////////////////////////////////////////////////////////////

void WriteInstantAndAvgSolutionHighOrder::setup()
{
CFAUTOTRACE;

  WriteSolution::setup();

}

//////////////////////////////////////////////////////////////////////////////

void WriteInstantAndAvgSolutionHighOrder::unsetup()
{
  CFAUTOTRACE;

  WriteSolution::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WriteInstantAndAvgSolutionHighOrder::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = WriteSolution::needsSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WriteInstantAndAvgSolutionHighOrder::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_pastAvgStates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace TecplotWriter

} // namespace COOLFluiD
