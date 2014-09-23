#ifndef COOLFluiD_Numerics_LMaestro_hh
#define COOLFluiD_Numerics_LMaestro_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/Maestro.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
namespace LoopMaestro {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a LMaestro. It controls the flow of actions in the SubSystem's through Event's.
/// @author Tiago Quintino
/// @author Thomas Wuilbaut
class LMaestro : public Framework::Maestro {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments.
  LMaestro(const std::string& name);

  /// Default destructor.
  ~LMaestro();

  /// Takes control of the simulation
  Common::Signal::return_t control ( Common::Signal::arg_t );

private: //data

  /// list of initial files
  std::vector<std::string> m_init_files;

}; // end of class LMaestro

//////////////////////////////////////////////////////////////////////////////

} // namespace LoopMaestro
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_LMaestro_hh
