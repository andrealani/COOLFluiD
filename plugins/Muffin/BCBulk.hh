#ifndef COOLFluiD_Muffin_BCBulk_hh
#define COOLFluiD_Muffin_BCBulk_hh

#include "Muffin/BC.hh"

namespace COOLFluiD {
  namespace Muffin {

class System;


/// Bulk concentrations boundary condition
class BCBulk : public BC {

 public:  // functions

  /// Boundary condition constructor
  BCBulk(const std::string& name);

  /// Boundary condition destructor
  ~BCBulk()
  {}

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Apply boundary condition
  void applyOnSystemMITReM(const Common::SafePtr< System > s, const Common::SafePtr< Framework::TopologicalRegionSet > t);


 private:  // data (user non-configurable)

  /// If m_concentrations has been adjusted to be electrically neutral
  bool m_concentrations_recalculated;


 private:  // data (user configurable)

  /// Factor to multiply bulk concentrations with
  double m_factor;

  /// Concentrations to apply (overrides "Factor")
  std::vector< double > m_concentrations;

};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_BCBulk_hh

