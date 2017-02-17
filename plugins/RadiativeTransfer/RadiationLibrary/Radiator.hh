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

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  static std::string getClassName() { return "Radiator"; }

  void configure(Config::ConfigArgs& args){;}

  /// Constructor without arguments
  Radiator(const std::string& name);
  ~Radiator(){;}
  
  /// @return the number of bins for spectra reduction
  CFuint getNbBins() const {return m_nbBins;}
  
  /// @return the number of bands for spetra reduction 
  CFuint getNbBands() const {return m_nbBands;}
  
  virtual void setup();
  
  virtual void unsetup() {}
  
  virtual void setupSpectra(CFreal wavMin, CFreal wavMax) = 0;

  virtual CFreal getEmission( CFreal lambda, RealVector &s_o ) = 0;

  virtual CFreal getAbsorption( CFreal lambda, RealVector &s_o ) = 0;

  virtual CFreal getSpectraLoopPower() = 0;

  virtual void computeEmissionCPD() = 0;

  virtual void getRandomEmission(CFreal &lambda, RealVector &s_o ) = 0;

  void setRadPhysicsPtr(RadiationPhysics *radPhysicsPtr) {
    m_radPhysicsPtr = radPhysicsPtr;
  }

  void setRadPhysicsHandlerPtr(RadiationPhysicsHandler *radPhysicsHandlerPtr){
    m_radPhysicsHandlerPtr = radPhysicsHandlerPtr;
  }
  
  /// get the volume of the current cell
  CFreal getCurrentCellVolume() const;
  
  /// get the area of the current wall face
  CFreal getCurrentWallArea() const;
  
  /// get the cell volume
  CFreal getCellVolume(CFuint stateID) const;
  
  /// get the wall face area
  CFreal getWallArea(CFuint wallGeoID) const;
  
protected:
  
  const CFreal m_angstrom; 
  RadiationPhysics *m_radPhysicsPtr;
  RadiationPhysicsHandler *m_radPhysicsHandlerPtr;
  RandomNumberGenerator m_rand;
  
  /// array of state vectors
  Framework::DataHandle<Framework::State*, Framework::GLOBAL> m_states;
  
  /// array of cell volumes
  Framework::DataHandle<CFreal> m_volumes;

  /// array fo face areas
  Framework::DataHandle<CFreal> m_faceAreas;
  
  /// array of face centers
  Framework::DataHandle<CFreal> m_faceCenters;
  
  /// number of bands
  CFuint m_nbBands;
  
  /// number of bins
  CFuint m_nbBins;
};
  
//////////////////////////////////////////////////////////////////////////////

}
}

#endif
