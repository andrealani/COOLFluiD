#ifndef COOLFluiD_Muffin_SystemTemp_hh
#define COOLFluiD_Muffin_SystemTemp_hh

#include "Muffin/System.hh"

namespace COOLFluiD {
  namespace Muffin {


/// Class extending System for convection-diffusion of a scalar
class SystemTemp : public System {

 public:  // functions

  /// Temperature system constructor
  SystemTemp(const std::string& name);

  /// Temperature system destructor
  ~SystemTemp();

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Configure the command
  void configure(Config::ConfigArgs& args);

  /// Set private data before processing phase
  void setup();

  /// Assemble elements residuals
  void executeOnTrs();

  /// Get diffusivity coefficients for variables system is responsible for
  std::vector< double > getDiffusivity() {
    return std::vector< double >(Nsys,m_diffusivity);
  }

 private:  // data

  /// Diffusivity constant
  double m_diffusivity;

};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_SystemTemp_hh

