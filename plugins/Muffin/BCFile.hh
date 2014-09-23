#ifndef COOLFluiD_Muffin_BCFile_hh
#define COOLFluiD_Muffin_BCFile_hh

#include "Muffin/BC.hh"

namespace COOLFluiD {
  namespace Muffin {

class System;


/// File boundary condition
class BCFile : public BC {

 public:  // functions

  /// Boundary condition constructor
  BCFile(const std::string& name);

  /// Boundary condition destructor
  ~BCFile();

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Configure the boundary condition
  void configure(Config::ConfigArgs& args);

  /// Set private data before processing phase
  void setup();

  /// Necessary instructions before destructor call
  void unsetup();

  /// Apply boundary condition (generic for all equation systems)
  void apply(const Common::SafePtr< System > s);


 public:  // data (user-configurable)

  /// File to read
  std::string m_file_str;

  /// Names of variables to apply to (default all of PhysicalModel)
  std::vector< std::string > m_applyeqs_str;


 public:  // data (non user-configurable)

  /// Indices of variables to apply to
  std::vector< int > m_applyeqs;

};


  }  // namespace Muffin
}  // namespace COOLFluiD

#endif // COOLFluiD_Muffin_BCFile_hh

