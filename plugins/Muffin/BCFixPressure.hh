#ifndef COOLFluiD_Muffin_BCFixPressure_hh
#define COOLFluiD_Muffin_BCFixPressure_hh

#include "Muffin/BC.hh"

namespace COOLFluiD {
  namespace Muffin {

class System;


/// Fixed pressure boundary condition
class BCFixPressure : public BC {

 public:  // functions

  /// System constructor
  BCFixPressure(const std::string& name);

  /// System destructor
  ~BCFixPressure();

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Set private data before processing phase
  void setup();

  /// Apply boundary condition
  void applyOnSystemFlow(const Common::SafePtr< System > s, const Common::SafePtr< Framework::TopologicalRegionSet > t);


 public:  // data

  /// Value
  double m_value;

};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_BCFixPressure_hh

