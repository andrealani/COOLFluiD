#ifndef COOLFluiD_Muffin_BCFunction_hh
#define COOLFluiD_Muffin_BCFunction_hh

#include "Muffin/BC.hh"

namespace COOLFluiD {
  namespace Muffin {

class System;


/// Function boundary condition
class BCFunction : public BC {

 public:  // functions

  /// Boundary condition constructor
  BCFunction(const std::string& name);

  /// Boundary condition destructor
  ~BCFunction();

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Configure the boundary condition
  void configure(Config::ConfigArgs& args);

  /// Apply boundary condition (generic for all equation systems)
  void apply(const Common::SafePtr< System > s);

  /// Set private data before processing phase
  void setup();


 private:  // functions

  /// Parse VectorialFunction
  void parse();


 private:  // data (user-configurable)

  /// Definition of the functions
  std::vector< std::string > m_function_str;

  /// Names of variables to apply to
  std::vector< std::string > m_applyvars_str;


 private:  // data (non user-configurable)

  /// VectorialFunction to use
  Framework::VectorialFunction m_function;

  /// Indices of variables to apply to
  std::vector< int > m_applyvars;
};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_BCFunction_hh

