#ifndef TemplateClassLoader_hh
#define TemplateClassLoader_hh

#ifdef CF_HAVE_CONFIG_H
  #include "coolfluid_config.h"
#endif

#if CF_HAVE_LIBTOOL
  #include <ltdl.h> // LibTool Dynamic Loading library header
#endif // CF_HAVE_LIBTOOL

#include <sstream>
#include <fstream>
#include <set>


#include "Common/CFLog.hh"
#include "Environment/DirPaths.hh"
#include "Environment/Factory.hh"
#include "Common/SelfRegistPtr.hh"
#include "AutoTemplateLoader/CodeBuilders.hh"

template < typename BASE, typename COMPILER = GCC >
class TemplateClassLoader {

  public: // typedefs

    /// the provider type for the BASE
    typedef typename BASE::PROVIDER ConcreteProvider;

    /// the type for the collection of parameters
    typedef typename std::vector< std::pair< std::string,std::string > > ParamsType;

  public: // functions

    TemplateClassLoader()
    {
      path = COOLFluiD::Environment::DirPaths::getInstance().getBaseDir().string();
    }

    COMPILER& Compiler() { return compiler; }

    COOLFluiD::Common::SafePtr<ConcreteProvider> loadClass (const std::string& classname, const ParamsType& params)
    {
      using namespace COOLFluiD;

      std::pair<std::string,std::string> name;

      try {
        name = create_lib(classname,params);
      }
      catch (std::string e) {
        CFLog(ERROR, "library creation exception : " << e << "\n");
      }

      std::string& libname  = name.first;
      std::string& concretename = name.second;

      try {
        load_module(libname);
      }
      catch (std::string e) {
        CFLog(ERROR, "module load exception : " << e << "\n");
      }

      COOLFluiD::Common::SafePtr<ConcreteProvider> prov = COOLFluiD::Environment::Factory<BASE>::getInstance().getProvider(concretename);
      m_loadedClasses[classname].insert(make_paramsname(params));
      return prov;
    }

    COOLFluiD::Common::SelfRegistPtr<BASE> createObject (const std::string& classname, const ParamsType& params)
    {
      std::string paramsname = make_paramsname(params);
      if (m_loadedClasses[classname].count(paramsname) != 0)
      {
        COOLFluiD::Common::SafePtr<ConcreteProvider> prov = COOLFluiD::Environment::Factory<BASE>::getInstance().getProvider(make_concretename(classname,paramsname));
        return prov->create();
      }
      else
      {
        std::string msg;
        msg = BASE::getClassName()
            + " was not yet instantiated with concrete class "
            + classname
            + " and parametrized with parameters "
            + paramsname;
        throw msg;
      }
    }

    private: // functions

      void load_module(const std::string& mod)
      {
        using namespace COOLFluiD;

        lt_dlsetsearchpath(path.c_str());
        lt_dlhandle hdl = lt_dlopenext(mod.c_str());
        if (hdl != 0)
        {
          CFLog(VERBOSE, "load_module: loaded " << mod << "\n");;
        }
        else
        {
          std::string error(lt_dlerror());
          throw error;
        }
      }

      std::string make_paramsname(const ParamsType& params)
      {
        std::string result;
        typename ParamsType::const_iterator itr = params.begin();
//         if (itr != params.end())
//         {
//           result += (*itr).first;
//           ++itr;
//         }
        for (; itr != params.end(); ++itr)
        {
          result += "_" + (*itr).first;
        }
        return result;
      }

      std::string make_templateparamsname(const ParamsType& params)
      {
        std::string result;
        typename ParamsType::const_iterator itr = params.begin();
        if (itr != params.end())
        {
          result += (*itr).first;
          ++itr;
        }
        for (; itr != params.end(); ++itr)
        {
          result += " , " + (*itr).first;
        }
        return result;
      }

      std::string make_concretename(const std::string& classname, const std::string& param)
      {
        return classname + param;
      }

      std::string make_libname(const std::string& classname, const std::string& param)
      {
        return compiler.libprefix() + make_concretename(classname,param);
      }

      std::string make_headerfilename(const std::string& classname)
      {
        return classname + ".hh";
      }

      std::string make_sourcefilename(const std::string& classname, const std::string& param)
      {
        return make_concretename(classname,param) + ".cxx";
      }

      std::string make_objectfilename(const std::string& classname, const std::string& param)
      {
        return make_concretename(classname,param) + compiler.objectext();
      }

      std::pair<std::string,std::string> create_lib (const std::string& classname, const ParamsType& params)
      {
        using namespace COOLFluiD;
        std::string paramsname = make_paramsname(params);

        CFLog(NOTICE, "CLASS  : " << classname  << " PARAMS : " << paramsname << "\n");

        std::string concretename = make_concretename(classname,paramsname);
        std::string libname  = make_libname(classname,paramsname);

        std::string headerfile = make_headerfilename(classname);
        std::string sourcefile = make_sourcefilename(classname,paramsname);

        std::string objfile = make_objectfilename(classname,paramsname);

        std::string baseclass = BASE::getClassName();

        std::ofstream fout(sourcefile.c_str());

        typename ParamsType::const_iterator itr = params.begin();
        for (; itr != params.end(); ++itr)
        {
          if (!(*itr).second.empty()) { fout << "#include \"" << (*itr).second << "\"" << std::endl; }
        }
        fout << "#include <cassert>" << std::endl;
        fout << "#include \"AutoTemplateLoader/AutoTemplateLoader.hh\"" << std::endl;
        fout << "#include \"" << headerfile << "\"" << std::endl;
        fout << "#include \"Environment/ObjectProvider.hh\"" << std::endl;
        fout
            << " COOLFluiD::Environment::ObjectProvider< "
            << classname << " < " << make_templateparamsname(params) << " > "
            << " , "
            << baseclass
            << " , "
            << "COOLFluiD::AutoTemplateLoader::AutoTemplateLoaderModule"
            << " > "
            << " a"<< classname << paramsname << "provider "
            << " ( \"" << concretename << "\" )"
            << " ; "
            << std::endl;

        fout.close();

        std::ostringstream ccom;
        ccom << compiler.compilecom()
            << sourcefile;
        system_execute(ccom.str());

        std::ostringstream lcom;
        lcom << compiler.linkcom()
            << libname << compiler.libext()
            << compiler.linkopts() << " "
            << objfile;
        system_execute(lcom.str());

        return std::make_pair<std::string,std::string>(libname,concretename);
      }

      int system_execute(const std::string& call)
      {
        using namespace COOLFluiD;
        CFLog(NOTICE, "Executing " << call << "\n");
#if 1
        return system(call.c_str());
#else
        return 0;
#endif
      }

      private: // data

        std::string path;
        std::map < std::string, std::set <std::string> > m_loadedClasses;
        COMPILER compiler;

}; // end TemplateClassLoader

#endif
