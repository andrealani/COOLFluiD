#ifndef COOLFluiD_RadiativeTransfer_GreyRadiator_hh
#define COOLFluiD_RadiativeTransfer_GreyRadiator_hh

#include "RadiativeTransfer/RadiationLibrary/Radiator.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

class GreyRadiator : public Radiator
{
public:

  static std::string getClassName() { return "GreyRadiator"; }
  static void defineConfigOptions(Config::OptionList& options);
  void configure(Config::ConfigArgs& args);

  /// Constructor without arguments
  GreyRadiator(const std::string& name);
  ~GreyRadiator(){;}

  void setup();

  void setupSpectra(CFreal wavMin, CFreal wavMax);

  CFreal getEmission( CFreal lambda, RealVector &s_o );

  CFreal getAbsorption( CFreal lambda, RealVector &s_o );

  CFreal getSpectraLoopPower();

  void computeEmissionCPD(){;}

  void getRandomEmission(CFreal &lambda, RealVector &s_o );
  
  void getData(){}
  
protected:
  
  virtual CFreal getCurrentElemTemperature() {return getCurrentCellTemperature();}
  
  virtual CFreal getSpaceIntegrator(){
    return getCurrentCellVolume() * 4.; //4 is from the solid angle integration
  }
  
  virtual void getRandomDirections(CFuint dim, RealVector &s_o ){
    getSphericalDirections(dim, s_o);
  }
  
  void getHemiDirections(CFuint dim, RealVector &s_o);
  
  void getSphericalDirections(CFuint dim, RealVector &s_o){
    m_rand.sphereDirections(dim, s_o);
  }

  CFreal computeComulativePlankFraction(CFreal lambda, CFreal T);
  inline CFreal computePlank(const CFreal lambda, const CFreal T);
  inline CFreal computeStefanBoltzmann(const CFreal T);
  CFreal getCurrentCellTemperature();
  CFreal getCurrentWallTemperature();

  CFreal getCellVolume();
  CFreal getFaceArea();

  CFreal m_minWav;
  CFreal m_maxWav;

  CFuint m_tempID;

  CFreal m_emsCoeff;
  CFreal m_absCoeff;

  bool m_allIsGrey;
  CFuint m_nbComPlankTerms;
  CFuint m_TRStypeID;
  static const CFreal m_Planck          = 6.6260693e-34 ;
  static const CFreal m_Boltzmann       = 1.3806580e-23 ;
  static const CFreal m_StefanBoltzmann = 5.67037321e-08;
  static const CFreal m_SpeedOfLight    = 2.99792458e8  ;
  static const CFreal m_Angstrom        = 1.0000000e-10 ;

  //m_c1=std::pow(m_Boltzmann,4)/(std::pow(m_Planck,3)*std::pow(m_SpeedOfLight,2))/m_StefanBoltzmann/(2*m_Pi);
  //m_c2=m_Planck*m_SpeedOfLight/m_Boltzmann;
  //m_c3=2.*m_Planck*std::pow(m_SpeedOfLight,2);
  static const CFreal m_c1 = 1.5399311364756e-01;
  static const CFreal m_c2 = 1.4387673140816e-02;
  static const CFreal m_c3 = 1.1910428196088e-16;
  
protected:
  
  /// array of face normals
  Framework::DataHandle<CFreal> m_socketNormals;
  
  /// problem dimension
  CFuint m_dim2;
  
  /// array storing temporary normal 
  RealVector m_normal;
};
  

/// /////////////////////////////////////////////////////////////////
/// \brief The GreyWallRadiator class
///
class GreyWallRadiator : public GreyRadiator
{
public:
  
  GreyWallRadiator(const std::string& name): GreyRadiator(name) {}
  
  ~GreyWallRadiator(){}
  
  static std::string getClassName() { return "GreyWallRadiator"; }

protected:
  virtual inline CFreal getCurrentElemTemperature(){
    return getCurrentWallTemperature();
  }

  virtual inline CFreal getSpaceIntegrator(){
    return getCurrentWallArea();
  }

  virtual inline void getRandomDirections(CFuint dim, RealVector &s_o ){
    getHemiDirections(dim, s_o);
  }
};

//////////////////////////////////////////////////////////////////////////////

}
}
#endif
