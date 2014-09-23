#include "RadiativeTransfer/RadiationLibrary/Reflector.hh"
namespace COOLFluiD {

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

class SpecularReflector : public Reflector
{
public:

  static std::string getClassName() { return "SpecularReflector"; }
  static void defineConfigOptions(Config::OptionList& options){}
  void configure(Config::ConfigArgs& args){}

  SpecularReflector(const std::string& name);
  ~SpecularReflector(){;}

  virtual void setup(){}

  void setupSpectra(CFreal wavMin, CFreal wavMax){}

  CFreal getReflectionProbability( CFreal lambda, RealVector &s_o, RealVector &s_i, RealVector &normal);

  void computeReflectionCPD(){}

  void getRandomDirection(CFreal &lambda, RealVector &s_o, RealVector &s_i, RealVector &normal );

  void getData(){}

private:
};

}
}
