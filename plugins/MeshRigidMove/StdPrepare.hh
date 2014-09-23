#ifndef COOLFluiD_Numerics_MeshRigidMove_StdPrepare_hh
#define COOLFluiD_Numerics_MeshRigidMove_StdPrepare_hh

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
class StdPrepare : public RigidMoveCom {
public:

  /**
   * Constructor.
   */
  explicit StdPrepare(std::string name) : RigidMoveCom(name)
  {
  }

  /**
   * Destructor.
   */
  ~StdPrepare()
  {
  }

  /**
   * Execute Processing actions
   */
  void execute();

private: // data

}; // class StdPrepare

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshRigidMove

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshRigidMove_StdPrepare_hh

