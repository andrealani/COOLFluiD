#include "Environment/FileHandlerInputConcrete.hh"
#include "Environment/ObjectProvider.hh"
#include "Environment/DirectFileAccess.hh"
#include "Environment/Environment.hh"
#include "Environment/FileHandlerOutputConcrete.hh"
#include "Environment/DirectFileWrite.hh"
#include "Framework/SMaestro.hh"
#include "Framework/Framework.hh"

namespace COOLFluiD {
  
  class PluginsRegister {
  public:
    
    void registerAll() 
    {
      using namespace Environment;
      using namespace Framework;
      
      /*  typedef Environment::ObjectProvider<Environment::FileHandlerInputConcrete< Environment::DirectFileAccess >, 
	  Environment::FileHandlerInput,
	  Environment::EnvironmentModule > FHI_PROVIDER; 
	  Factory<FileHandlerInput>::getInstance().regist(new FHI_PROVIDER("DirectFileAccess"));
	  
	  typedef Environment::ObjectProvider<Environment::FileHandlerOutputConcrete< DirectFileWrite >, 
	  Environment::FileHandlerOutput, 
	  Environment::EnvironmentModule > FHO_PROVIDER; 
	  Factory<FileHandlerOutput>::getInstance().regist(new FHO_PROVIDER("DirectFileWrite")); 
	  
	  Factory<Maestro>::getInstance().regist(new  
	  ObjectProvider<SMaestro, Maestro, FrameworkLib, 1>("SimpleMaestro")); */
      
      Factory<FileHandlerInput>::getInstance().regist
	(new Environment::ObjectProvider<Environment::FileHandlerInputConcrete< Environment::DirectFileAccess >, 
	 Environment::FileHandlerInput,
	 Environment::EnvironmentModule >("DirectFileAccess"));
      
      Factory<FileHandlerOutput>::getInstance().regist
	(new Environment::ObjectProvider<Environment::FileHandlerOutputConcrete< DirectFileWrite >, 
	 Environment::FileHandlerOutput, 
	 Environment::EnvironmentModule >("DirectFileWrite")); 
      
      Factory<Maestro>::getInstance().regist
	(new ObjectProvider<SMaestro, Maestro, FrameworkLib, 1>("SimpleMaestro"));
    }
  };
}
