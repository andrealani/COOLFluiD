#ifndef COOLFluiD_RadiativeTransfer_NullRadiator_hh
#define COOLFluiD_RadiativeTransfer_NullRadiator_hh

#include "RadiativeTransfer/RadiationLibrary/Radiator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

class NullRadiator : public Radiator
{
public:

  static std::string getClassName() { return "NullRadiator"; }
  static void defineConfigOptions(Config::OptionList& options){;}
  void configure(Config::ConfigArgs& args){;}

  /// Constructor without arguments
  NullRadiator(const std::string& name):
    Radiator(name)
  {
    addConfigOptionsTo(this);
  }

  ~NullRadiator(){;}

  void setup(){;}

  void setupSpectra(CFreal wavMin, CFreal wavMax){;}

  CFreal getEmission( CFreal lambda, RealVector &s_o ){ return 0.; }

  CFreal getAbsorption( CFreal lambda, RealVector &s_o ){ return 0.; }

  CFreal getSpectaLoopPower(){return 0.;}

  void computeEmissionCPD(){;}

  void getRandomEmission(CFreal &lambda, RealVector &s_o ){
    CFLog(INFO,"Called Emission form NullRadiator");
  }

  void getData(){;}
};

}
}

#endif
