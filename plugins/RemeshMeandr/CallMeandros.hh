#ifndef COOLFluiD_Numerics_RemeshMeandros_CallMeandros_hh
#define COOLFluiD_Numerics_RemeshMeandros_CallMeandros_hh

//////////////////////////////////////////////////////////////////////////////

#include "RMeshMeData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RemeshMeandros {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a NumericalCommand action to be
/// sent to Domain to be executed in order to ComputeSpaceResidual the MeshData.
/// @author Jurek Majewski
class CallMeandros : public RMeshMeCom
{
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  explicit CallMeandros(const std::string& name);

  /// Destructor.
  ~CallMeandros();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  void setup();

  /// Execute Processing actions
  void execute();

private: // data

  /// reference lenght in mesh
  CFreal _refH;

}; // class CallMeandros

//////////////////////////////////////////////////////////////////////////////

    } // namespace RemeshMeandros

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_RemeshMeandros_CallMeandros_hh

