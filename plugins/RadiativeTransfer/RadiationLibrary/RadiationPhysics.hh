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

  /// Constructor
  RadiationPhysics(const std::string& name);
  
  /// Destructor
  ~RadiationPhysics(){}
  
  /// setup the spectra
  void setupSpectra(CFreal wavMin,CFreal wavMax);
  
  /// setup private data
  void setup();
  
  /// configure this object from given options
  void configure(Config::ConfigArgs& args);
  
  /// define the configurable options statically
  static void defineConfigOptions(Config::OptionList& options);
  
  /// @return the class name
  static std::string getClassName() { return "RadiationPhysics"; }
  
  /// @return the name of the TRS to which RadiationPhysics is applied
  std::string getTRSname(){return m_TRSname;}
  
  /// @return the type of the TRS to which RadiationPhysics is applied
  std::string getTRStype(){return m_TRStype;}
  
  /// @return the type ID of the TRS to which RadiationPhysics is applied
  TRStypeID getTRStypeID(){return m_TRStypeID;}
  
  /// @return the radiation physics itself
  void getRadPhysicsPtr(RadiationPhysics* ptr) { ptr = this; }
  
  /// get the local IDs of the Wall TRS ghost states and the face
  /// @param statesID    local IDs of the Wall TRS ghost states
  /// @param wallGeoIdx  local IDs of the Wall TRS faces
  void getWallStateIDs(std::vector<CFuint> &statesID,
                       std::vector<CFuint>& wallGeoIdx);
  
  /// @return the cell state ID
  void getCellStateIDs(std::vector<CFuint> &statesID);
  
  /// @return the wall state
  Common::SafePtr<Framework::State> getWallState(CFuint wallTrsIdx){return &(m_interpolatedStates[wallTrsIdx]);}
  
  /// set the @see RadiationPhysicsHandler
  void setRadPhysicsHandlerPtr( RadiationPhysicsHandler* ptr ) {m_radPhysicsHandlerPtr = ptr;}
  
  /// @return the pointer to the radiator
  Common::SafePtr<Radiator> getRadiatorPtr() const {return m_radiator.getPtr();}
  
  /// @return the pointer to the reflector
  Common::SafePtr<Reflector> getReflectorPtr() const {return m_reflector.getPtr();}
  
  /// compute the interpolated states
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
