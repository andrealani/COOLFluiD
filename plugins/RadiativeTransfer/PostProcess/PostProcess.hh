#ifndef COOLFluiD_Post_Process_hh
#define COOLFluiD_Post_Process_hh

#include <vector>

#include "FiniteVolume/CellCenterFVM.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/SocketBundleSetter.hh"

namespace COOLFluiD {

namespace RadiativeTransfer {


using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;


class PostProcess  : public Common::OwnedObject,
                     public Config::ConfigObject,
                     public Common::NonCopyable<PostProcess>,
                     public Framework::SocketBundleSetter {
public:

  typedef Environment::ConcreteProvider<PostProcess,1> PROVIDER;
  typedef const std::string& ARG1;

  /**
  * Defines the Config Option's of this class
  * @param options a OptionList where to add the Option's
  */
  static void defineConfigOptions(Config::OptionList& options);

  PostProcess(const string &name);

  void configure ( Config::ConfigArgs& args );

  /// Default destructor
  virtual ~PostProcess();

  /// Gets the Class name
  static std::string getClassName() { return "PostProcess"; }

  virtual void runPostProcess(DataHandle<CFreal> dataVector) = 0;

}; // end of class


//////////////////////////////////////////////////////////////////////////////

} // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Post_Process_hh
