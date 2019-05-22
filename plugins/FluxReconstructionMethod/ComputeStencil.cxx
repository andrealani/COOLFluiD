#include "FluxReconstructionMethod/ComputeStencil.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

void ComputeStencil::defineConfigOptions(Config::OptionList& options)
{ 
  options.addConfigOption< std::vector<std::string> >
    ("TRSNames","Names of the TRSs to which the stencil computation is applied");
}
  
//////////////////////////////////////////////////////////////////////////////
      
ComputeStencil::ComputeStencil(const std::string& name) :
  Common::OwnedObject(),
  ConfigObject(name),
  Common::NonCopyable<ComputeStencil>(),
  socket_states("Null"),
  socket_nodes("Null"),
  socket_stencil("Null")
{
  addConfigOptionsTo(this);
  
  _trsNames = std::vector<std::string>();
  setParameter("TRSNames",&_trsNames);
}

//////////////////////////////////////////////////////////////////////////////
  
ComputeStencil::~ComputeStencil()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeStencil::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
