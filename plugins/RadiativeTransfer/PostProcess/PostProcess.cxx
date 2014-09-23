#include "PostProcess.hh"

namespace COOLFluiD {

namespace RadiativeTransfer {

PostProcess::PostProcess(const std::string& name) :
  Common::OwnedObject(),
  ConfigObject(name),
  Common::NonCopyable<PostProcess>(),
  Framework::SocketBundleSetter()
{
  addConfigOptionsTo(this);
}

PostProcess::~PostProcess()
{
}

void PostProcess::defineConfigOptions(Config::OptionList& options)
{
}

void PostProcess::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);
}




}

}
