#ifndef COOLFluiD_Numerics_TAdptMaestro_hh
#define COOLFluiD_Numerics_TAdptMaestro_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/Maestro.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RemeshMeandros {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a TAdptMaestro. It controls the flow of actions in the SubSystem's through Event's.
/// @author Jurek Majewski
class TAdptMaestro : public Framework::Maestro {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor without arguments.
  TAdptMaestro(const std::string& name);

  /// Default destructor.
  ~TAdptMaestro();

  /// Takes control of the simulation
  Common::Signal::return_t control ( Common::Signal::arg_t );

private: //data

  std::vector<std::string> _initialFilesStr;

}; // end of class TAdptMaestro

//////////////////////////////////////////////////////////////////////////////

    } // namespace RemeshMeandros

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_TAdptMaestro_hh
