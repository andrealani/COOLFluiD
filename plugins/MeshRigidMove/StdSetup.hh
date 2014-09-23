#ifndef COOLFluiD_Numerics_MeshRigidMove_StdSetup_hh
#define COOLFluiD_Numerics_MeshRigidMove_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "RigidMoveData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshRigidMove {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to Domain to be executed in order to setup the MeshData.
   *
   * @author Thomas Wuilbaut
   *
   */
class StdSetup : public RigidMoveCom {
public:

  /**
   * Constructor.
   */
  explicit StdSetup(std::string name) : RigidMoveCom(name)
  {
  }

  /**
   * Destructor.
   */
  ~StdSetup()
  {
  }

  /**
   * Execute Processing actions
   */
  void execute();


}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshRigidMove

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshRigidMove_StdSetup_hh

