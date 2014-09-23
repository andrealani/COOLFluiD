#ifndef COOLFluiD_Numerics_MeshFEMMove_StdPrepare_hh
#define COOLFluiD_Numerics_MeshFEMMove_StdPrepare_hh

//////////////////////////////////////////////////////////////////////////////

#include "FEMMoveData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshFEMMove {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to Domain to be executed in order to setup the MeshData.
   *
   * @author Thomas Wuilbaut
   *
   */
class StdPrepare : public FEMMoveCom {
public:

  /**
   * Constructor.
   */
  explicit StdPrepare(std::string name) : FEMMoveCom(name)
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

    } // namespace MeshFEMMove

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshFEMMove_StdPrepare_hh

