#ifndef CodeBuilders_hh
#define CodeBuilders_hh



class BaseCompiler
{
  public:

    BaseCompiler() : includeflags(), linkerflags()
    {
    };

  virtual std::string libprefix  () = 0;
  virtual std::string libext     () = 0;
  virtual std::string objectext  () = 0;
  virtual std::string compilecom () = 0;
  virtual std::string linkcom    () = 0;
  virtual std::string linkopts   () = 0;

  void addToIncludeFlags(std::string flags) { includeflags += flags; };
  void addToLinkerFlags (std::string flags) { linkerflags  += flags; };

  protected: // data
    std::string includeflags;
    std::string linkerflags;
};

#if CF_HAVE_LIBTOOL
struct GCCLibTool : public BaseCompiler
{
  std::string libprefix  () { return "lib"; }
  std::string libext     () { return ".la"; }
  std::string objectext  () { return ".lo"; }
  std::string compilecom () { return std::string("libtool --tag=CXX --mode=compile g++ ") + includeflags + std::string(" -c "); }
  std::string linkcom    () { return std::string("libtool --tag=CXX --mode=link g++ ") + linkerflags + std::string(" -shared -o "); }
  std::string linkopts   () { return " -rpath `pwd` -no-undefined -module -avoid-version "; }
};
#endif // CF_HAVE_LIBTOOL

struct GCC : public BaseCompiler
{
  std::string libprefix  () { return "lib"; }
  std::string libext     () { return ".so"; }
  std::string objectext  () { return ".o"; }
  std::string compilecom () { return std::string("g++ -fPIC ") + includeflags + std::string(" -c "); }
  std::string linkcom    () { return std::string("g++ -fPIC ") + linkerflags + std::string(" -shared -o "); }
  std::string linkopts   () { return ""; }
};

struct LAMGCC : public BaseCompiler
{
  std::string libprefix  () { return "lib"; }
  std::string libext     () { return ".so"; }
  std::string objectext  () { return ".o"; }
  std::string compilecom () { return std::string("mpic++ -fPIC ") + includeflags + std::string(" -c "); }
  std::string linkcom    () { return std::string("mpic++ -fPIC ") + linkerflags + std::string(" -shared -o "); }
  std::string linkopts   () { return ""; }
};

#endif
