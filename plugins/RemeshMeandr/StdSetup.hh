#ifndef COOLFluiD_Numerics_RemeshMeandros_StdSetup_hh
#define COOLFluiD_Numerics_RemeshMeandros_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "RMeshMeData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RemeshMeandros {

//////////////////////////////////////////////////////////////////////////////

  /// This class represents a NumericalCommand action to be
  /// sent to Domain to be executed in order to setup the MeshData.
class StdSetup : public RMeshMeCom {
public:

  /// Constructor.
  explicit StdSetup(std::string name) : RMeshMeCom(name)
  {
  }

  /// Destructor.
  ~StdSetup()
  {
  }

  /// Execute Processing actions
  void execute();

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace RemeshMeandros

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_RemeshMeandros_StdSetup_hh

