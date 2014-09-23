#ifndef COOLFluiD_Catalycity_BasicCatalycityModel_hh
#define COOLFluiD_Catalycity_BasicCatalycityModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/CatalycityModel.hh"
#include "Framework/PhysicalConsts.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Catalycity {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a basic catalycity library
/// @author Andrea Lani
class BasicCatalycityModel : public Framework::CatalycityModel {
public:
  
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);
  
  /// Constructor without arguments
  BasicCatalycityModel(const std::string& name);
  
  /// Default destructor
  ~BasicCatalycityModel();
  
  /// set up private data
  void setup();
  
  /// Compute the mass production
  /// @param temp  roto-translational temperature
  /// @param rho   density [kg/m^3]
  /// @param ys    array of species mass fractions
  /// @param mp    array of species mass production due to catalycity [kg m^-3 s^-1]
  void computeMassProduction(const CFreal temp, 
			     const CFreal rho,
			     const RealVector& ys, 
			     RealVector& mp);
  
private:
  
  /// @author: Jan Thoemel
  /// @date: October 2006
  /// remark: computes mole production rate
  /// @param concentrationreactant-concentration of reactant in [mol/m^3]
  /// @param molarmassreactant-molecular mass,[kg/particle]
  /// @param Twall-temperature of the wall,[K]
  /// @param catalycity-recombination coefficient,[1]
  /// output: reactantconsumption-moleproduction at the wall, [mol/m^2]
  void catalycityMoleProduction (CFreal molarFractionReactant,
				 CFreal molarMassReactant,
				 CFreal Twall,
				 CFreal catalycity,
				 CFreal& reactantConsumption)
  {
    using namespace COOLFluiD::Framework;
    using namespace COOLFluiD::MathTools;
    
    const CFreal cat = 2.*catalycity/(2.-catalycity);
    reactantConsumption = -cat*molarFractionReactant*
      std::sqrt(PhysicalConsts::UnivRgas()*Twall/(2.*MathConsts::CFrealPi()*molarMassReactant));
  }
  
private:
  
  /// species molar masses
  RealVector m_mmasses;
  
  /// species molar concentrations [mol/m^3] 
  RealVector m_mconcentrations;
  
  ///ID identifying the basic model to use
  CFuint m_modelID;
  
  ///catalycity coefficient per species
  std::vector<CFreal> m_catalycity;
  
}; // end of class BasicCatalycityModel
    
//////////////////////////////////////////////////////////////////////////////

    } // namespace Catalycity

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Catalycity_BasicCatalycityModel_hh
