#ifndef COOLFluiD_Muffin_BCFixVelocity_hh
#define COOLFluiD_Muffin_BCFixVelocity_hh

#include "Muffin/BC.hh"

namespace COOLFluiD {
  namespace Muffin {

class System;


/// Fixed velocity boundary condition
class BCFixVelocity : public BC {

public:  // functions

  /// System constructor
  BCFixVelocity(const std::string& name);

  /// System destructor
  ~BCFixVelocity();

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Configure the command
  void configure(Config::ConfigArgs& args);

  /// Set private data before processing phase
  void setup();

  /// Necessary instructions before destructor call
  void unsetup();

  /// Apply boundary condition
  void applyOnSystemFlow(const Common::SafePtr< System > s, const Common::SafePtr< Framework::TopologicalRegionSet > t);
  void applyOnSystemTemp(const Common::SafePtr< System > s, const Common::SafePtr< Framework::TopologicalRegionSet > t);
  void applyOnSystemTurb(const Common::SafePtr< System > s, const Common::SafePtr< Framework::TopologicalRegionSet > t);


private:  // functions

  /// Read velocity profile from file and interpolate
  void readFile(const std::string& filename);


public:  // data

  /// Values vector
  std::vector< double > invals;

  /// Fixed velocity type
  std::string type_str;

  /// File to read if type is set as such
  std::string file_str;

};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_BCFixVelocity_hh

