#include "RemeshMeandr/RemeshMeandr.hh"
#include "CallMeandros.hh"
#include "Framework/MeshData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Environment/DirPaths.hh"
#include "Common/OSystem.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RemeshMeandros {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<CallMeandros, RMeshMeData, RemeshMeandrModule> CallMeandrosProvider("StdCallMeandros");

//////////////////////////////////////////////////////////////////////////////

void CallMeandros::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("RefH","Reference H.");
}

//////////////////////////////////////////////////////////////////////////////

CallMeandros::CallMeandros(const std::string& name) :
RMeshMeCom(name)
{
   addConfigOptionsTo(this);
  _refH = 10.0;
   setParameter("RefH",&_refH);



  CFLog( ERROR, "inside CallMeandros constructor\n");
}

//////////////////////////////////////////////////////////////////////////////

CallMeandros::~CallMeandros()
{
}

//////////////////////////////////////////////////////////////////////////////

void CallMeandros::setup()
{
  CFLogDebugMin( "CallMeandros::setup() BEGIN\n");

  // first call parent method
  RMeshMeCom::setup();

  CFLogDebugMin( "CallMeandros::setup() END\n");
}

//////////////////////////////////////////////////////////////////////////////

void CallMeandros::execute()
{
  CFLog( ERROR, "CallMeandros::executeOnTrs() BEGIN\n" );

	CFout << " call meandros\n";

  std::string command;
  std::string resultsDir = Environment::DirPaths::getInstance().getResultsDir().string();
  std::string workingDir = Environment::DirPaths::getInstance().getWorkingDir().string();

  std::string meDir = Environment::DirPaths::getInstance().getWorkingDir().string() + getMethodData().getMeandrosDir();

  CFout << resultsDir << "\n";
  CFout << workingDir << "\n";
  CFout << meDir      << "\n";

  command = "cd " + meDir + "; ./generate.scr";
  CFout << command << "\n";

  Common::OSystem::getInstance().executeCommand( command);

  CFLogDebugMin( "CallMeandros::executeOnTrs() END\n");
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RemeshMeandros

  } // namespace Numerics

} // namespace COOLFluiD
