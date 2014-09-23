#ifndef BaseRunLib_hh
#define BaseRunLib_hh

#include "Common/OwnedObject.hh"
#include "Common/Stopwatch.hh"
#include "Environment/ConcreteProvider.hh"

using namespace COOLFluiD;

class BaseRunLib : public COOLFluiD::Common::OwnedObject
{
  public:

  typedef COOLFluiD::Environment::ConcreteProvider<BaseRunLib> PROVIDER;

   static std::string getClassName () { return "BaseRunLib"; }

  double test ( CFuint in_nbelems )
  {
      nbelems = in_nbelems;
      init();

      stp.start();
        compute();
      stp.stop();

      finalize();
      return stp.read();
  }

  double average_run ( CFuint in_nbelems, CFuint ntests )
  {
    double time = 0.;
    for ( CFuint i = 0; i < ntests; ++i)
    {
      time += test(in_nbelems);
    }
    return time / ntests;
  }

  protected: // hook methods

    virtual void compute  ( ) = 0;
    virtual void init     ( ) = 0;
    virtual void finalize ( ) = 0;

   protected: // data

    COOLFluiD::Common::Stopwatch<COOLFluiD::Common::WallTime> stp;
    unsigned int nbelems;

};

#endif // BaseRunLib_hh
