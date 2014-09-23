#include "RadiativeTransfer/RadiationLibrary/Reflector.hh"

namespace COOLFluiD {

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

class NullReflector: public Reflector
{
public:

  static std::string getClassName() { return "NullReflector"; }
  static void defineConfigOptions(Config::OptionList& options){}
  void configure(Config::ConfigArgs& args){}

  NullReflector(const std::string& name);
  ~NullReflector(){;}

  virtual void setup(){}

  void setupSpectra(CFreal wavMin, CFreal wavMax){}

  CFreal getReflectionProbability( CFreal lambda, RealVector &s_o, RealVector &s_i, RealVector &normal)
         { return 0.; }

  void computeReflectionCPD(){}

  void getRandomDirection(CFreal &lambda, RealVector &s_o, RealVector &s_i, RealVector &normal ){}

  void getData(){}

private:
};

}
}
