#include "Config/BiggerThanZero.hh"
#include "Maxwell/ConvMaxwellTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Config;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Maxwell {
      
//////////////////////////////////////////////////////////////////////////////

void ConvMaxwellTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal , Config::DynamicOption<> >("divBCleaningConst","Hyperbolic divB cleaning Constant");
  options.addConfigOption< CFreal , Config::DynamicOption<> >("divECleaningConst","Hyperbolic divE cleaning Constant");   
  options.addConfigOption< CFreal , Config::DynamicOption<> >("divBAdimCleaningConst","Adimensional hyperbolic divB cleaning Constant");       
  options.addConfigOption< CFreal > ("lightSpeedMax", "Speed of light. It can be reduced if it is still bigger than the speed of the fluid");
  //    options.addConfigOption< std::string >("correctionType","Name of correction for the projection scheme.");  
 
}         

//////////////////////////////////////////////////////////////////////////////

ConvMaxwellTerm::ConvMaxwellTerm(const std::string& name)
  : BaseTerm(name)
{
  addConfigOptionsTo(this);
  _divBCleaningConst = 1.;
  setParameter("divBCleaningConst",&_divBCleaningConst);    
  
  _divECleaningConst = 0.;
  setParameter("divECleaningConst",&_divECleaningConst); 
  
  _divBAdimCleaningConst = 1.;
  setParameter("divBAdimCleaningConst",&_divBAdimCleaningConst);
  
  _LightSpeed = 299792458.;		//light Speed m/s
  setParameter("lightSpeedMax",&_LightSpeed);
  
//   _correctionType = "Hyperbolic";
//   setParameter("correctionType",&_correctionType);   
   
   
//  _saveRate = 1;
//   setParameter("divBErrorFileSaveRate",&_saveRate);
  
//  _nameOutputFile = "divB";
//   setParameter("divBErrorFileName",&_nameOutputFile);
   
 
  //_electronMass = 5.4858e-7;		// Electron's mass [kg/mol] source:Standart Handbook for Electrical Engineerings 
  _electronMass = 9.1094e-31;		// Electron's mass [kg] source:Standart Handbook for Electrical Engineerings   
  _electronCharge = 1.602e-19;		// Electron's Charge [C]
  _protonMass = 1.6726e-27;		// Proton's mass [kg] source:Standart Handbook for Electrical Engineerings 
  //_protonMass = 1.00728e-3;		// Proton's mass [kg/mol] source:Standart Handbook for Electrical Engineerings
  _neutralMass = 1.6726e-27;		// Neutral's mass [kg/mol] equal to proton's mass   
  //_neutralMass = 1.00797e-3;		// Neutral's mass [kg/mol] source:Standart Handbook for Electrical Engineerings 
  _solarGravity = 274;			// Solar gravity [m/s^2] 
  
}
      
//////////////////////////////////////////////////////////////////////////////

ConvMaxwellTerm::~ConvMaxwellTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void ConvMaxwellTerm::configure ( Config::ConfigArgs& args )
{ 
  BaseTerm::configure(args);  
}

//////////////////////////////////////////////////////////////////////////////

void ConvMaxwellTerm::resizePhysicalData(RealVector& physicalData)
{
  // resize the physical data
  cf_assert(getDataSize() > 0);
  physicalData.resize(getDataSize());
}

//////////////////////////////////////////////////////////////////////////////

void ConvMaxwellTerm::setupPhysicalData()
{
  cf_assert(getDataSize() > 0);

  // set the size of each physical data in the StatesData
  resizePhysicalData(m_physicalData);
  resizePhysicalData(m_refPhysicalData);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace Maxwell

} // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
