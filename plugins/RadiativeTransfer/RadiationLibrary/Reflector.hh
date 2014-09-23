#ifndef COOLFluiD_RadiativeTransfer_Reflector_hh
#define COOLFluiD_RadiativeTransfer_Reflector_hh

#include "Common/OwnedObject.hh"
#include "Common/SetupObject.hh"
#include "Common/NonCopyable.hh"
#include "Environment/ConcreteProvider.hh"
#include "RadiativeTransfer/RadiativeTransferModule.hh"
#include "Framework/SocketBundleSetter.hh"
#include "RadiativeTransfer/Solvers/MonteCarlo/RandomNumberGenerator.hh"
#include "RadiationPhysics.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

class RadiationPhysicsHandler;
class RadiativePhysics;

class Reflector : public Common::OwnedObject,
                  public Config::ConfigObject
{
public:

  typedef Environment::ConcreteProvider<Reflector,1> PROVIDER;
  typedef const std::string& ARG1;

  static std::string getClassName() { return "Reflector"; }
  static void defineConfigOptions(Config::OptionList& options){;}
  void configure(Config::ConfigArgs& args){;}

  /// Constructor without arguments
  Reflector(const std::string& name);
  ~Reflector(){;}

  virtual void setup() = 0;

  virtual void setupSpectra(CFreal wavMin, CFreal wavMax) = 0;

  virtual CFreal getReflectionProbability( CFreal lambda, RealVector &s_o, RealVector &s_i, RealVector &normal ) = 0;

  virtual void computeReflectionCPD() = 0;

  virtual void getRandomDirection(CFreal &lambda, RealVector &s_o, RealVector &s_i, RealVector &normal ) = 0;

  virtual void getData() = 0;

  void setRadPhysicsPtr(RadiationPhysics *radPhysicsPtr) {
    m_radPhysicsPtr = radPhysicsPtr;
  }

  void setRadPhysicsHandlerPtr(RadiationPhysicsHandler *radPhysicsHandlerPtr){
    m_radPhysicsHandlerPtr = radPhysicsHandlerPtr;
  }

protected:
  RadiationPhysics *m_radPhysicsPtr;
  RadiationPhysicsHandler *m_radPhysicsHandlerPtr;
  RandomNumberGenerator m_rand;

};


}
}

#endif
