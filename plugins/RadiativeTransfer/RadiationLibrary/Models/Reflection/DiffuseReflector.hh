#include "RadiativeTransfer/RadiationLibrary/Reflector.hh"

namespace COOLFluiD {

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

class DiffuseReflector: public Reflector
{
public:

  static std::string getClassName() { return "DiffuseReflector"; }
  static void defineConfigOptions(Config::OptionList& options){}
  void configure(Config::ConfigArgs& args){}

  DiffuseReflector(const std::string& name);
  ~DiffuseReflector(){;}
  
  virtual void setup(){Reflector::setup();}
  
  void setupSpectra(CFreal wavMin, CFreal wavMax){}

  CFreal getReflectionProbability( CFreal lambda, RealVector &s_o, RealVector &s_i, RealVector &normal);

  void computeReflectionCPD(){}

  void getRandomDirection(CFreal &lambda, RealVector &s_o, RealVector &s_i, RealVector &normal );

};

}
}
