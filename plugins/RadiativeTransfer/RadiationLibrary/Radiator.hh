#ifndef COOLFluiD_RadiativeTransfer_Radiator_hh
#define COOLFluiD_RadiativeTransfer_Radiator_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/OwnedObject.hh"
#include "Common/SetupObject.hh"
#include "Common/NonCopyable.hh"
#include "Environment/ConcreteProvider.hh"
#include "RadiativeTransfer/RadiativeTransfer.hh"
#include "Framework/SocketBundleSetter.hh"
#include "RadiativeTransfer/Solvers/MonteCarlo/RandomNumberGenerator.hh"
#include "RadiativeTransfer/RadiationLibrary/RadiationPhysics.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

class RadiationPhysicsHandler;
class RadiationPhysics;

class Radiator : public Common::OwnedObject,
                 public Config::ConfigObject
{
public:

  typedef Environment::ConcreteProvider<Radiator,1> PROVIDER;
  typedef const std::string& ARG1;

  static std::string getClassName() { return "Radiator"; }
  static void defineConfigOptions(Config::OptionList& options){;}
  void configure(Config::ConfigArgs& args){;}

  /// Constructor without arguments
  Radiator(const std::string& name);
  ~Radiator(){;}

  virtual void setup() = 0;

  virtual void setupSpectra(CFreal wavMin, CFreal wavMax) = 0;

  virtual CFreal getEmission( CFreal lambda, RealVector &s_o ) = 0;

  virtual CFreal getAbsorption( CFreal lambda, RealVector &s_o ) = 0;

  virtual CFreal getSpectraLoopPower() = 0;

  virtual void computeEmissionCPD() = 0;

  virtual void getRandomEmission(CFreal &lambda, RealVector &s_o ) = 0;

  virtual void getData() = 0;

  void setRadPhysicsPtr(RadiationPhysics *radPhysicsPtr) {
    m_radPhysicsPtr = radPhysicsPtr;
  }

  void setRadPhysicsHandlerPtr(RadiationPhysicsHandler *radPhysicsHandlerPtr){
    m_radPhysicsHandlerPtr = radPhysicsHandlerPtr;
  }

  CFreal getCurrentCellVolume();
  CFreal getCurrentWallArea();

  CFreal getCellVolume( CFuint stateID );
  CFreal getWallArea(CFuint wallGeoID );

protected:
  
  const CFreal m_angstrom; 
  RadiationPhysics *m_radPhysicsPtr;
  RadiationPhysicsHandler *m_radPhysicsHandlerPtr;
  RandomNumberGenerator m_rand;
};
  
//////////////////////////////////////////////////////////////////////////////

}
}

#endif
