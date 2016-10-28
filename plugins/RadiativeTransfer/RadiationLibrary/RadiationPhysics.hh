#ifndef COOLFluiD_RadiativeTransfer_RadiationPhysics_hh
#define COOLFluiD_RadiativeTransfer_RadiationPhysics_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/OwnedObject.hh"
#include "Common/SetupObject.hh"
#include "Common/NonCopyable.hh"
#include "Environment/ConcreteProvider.hh"
#include "Environment/ObjectProvider.hh"
#include "RadiativeTransfer/RadiationLibrary/Radiator.hh"
#include "RadiativeTransfer/RadiationLibrary/RadiationPhysicsHandler.hh"
#include "RadiativeTransfer/RadiationLibrary/Reflector.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

class Radiator;
class Reflector;

namespace RadiativeDistTypes{
  enum TRStypeID {WALL=0, MEDIUM=1, BOUNDARY=2 };
  enum DistPhysicsID {ABSORPTION, EMISSION, SCATTERING, REFLECTION};
}

using namespace RadiativeDistTypes;

//////////////////////////////////////////////////////////////////////////////

class RadiationPhysics : public Common::OwnedObject,
                         public Config::ConfigObject 
{
public:

  typedef Environment::ConcreteProvider<RadiationPhysics,1> PROVIDER;
  typedef const std::string& ARG1;

  RadiationPhysics(const std::string& name);
  
  ~RadiationPhysics(){}
  
  void setupSpectra(CFreal wavMin,CFreal wavMax);
  
  void setup();
  
  void configure(Config::ConfigArgs& args);
  
  static void defineConfigOptions(Config::OptionList& options);
  
  static std::string getClassName() { return "RadiationPhysics"; }

  std::string getTRSname(){return m_TRSname;}
  
  std::string getTRStype(){return m_TRStype;}
  
  TRStypeID getTRStypeID(){return m_TRStypeID;}
  
  void getRadPhysicsPtr(RadiationPhysics* ptr) { ptr = this; }
  
  void getWallStateIDs(std::vector<CFuint> &statesID,
                       std::vector<CFuint>& wallGeoIdx);
  
  void getCellStateIDs(std::vector<CFuint> &statesID);
  
  Framework::State* getWallState(CFuint wallTrsIdx){return &(m_interpolatedStates[wallTrsIdx]);}
  
  void setRadPhysicsHandlerPtr( RadiationPhysicsHandler* ptr ) {m_radPhysicsHandlerPtr = ptr;}
  
  Radiator* const getRadiatorPtr(){return m_radiator.getPtr();}
  
  Reflector* const getReflectorPtr(){return m_reflector.getPtr();}
  
  void computeInterpolatedStates();

private:
  Common::SelfRegistPtr< Radiator > m_radiator;
  Common::SelfRegistPtr< Reflector > m_reflector;
  RadiationPhysicsHandler* m_radPhysicsHandlerPtr;
  TRStypeID m_TRStypeID;
  std::string m_radiatorName;
  std::string m_reflectorName;
  std::string m_TRSname;
  std::string m_TRStype;
  std::vector<Framework::State> m_interpolatedStates;
  std::vector<CFreal> m_faceAreas;
};

//////////////////////////////////////////////////////////////////////////////

} // end namespace RadiativeTransfer

}// end namespace COOLFluiD

#endif
