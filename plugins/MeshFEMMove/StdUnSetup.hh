#ifndef COOLFluiD_Numerics_MeshFEMMove_StdUnSetup_hh
#define COOLFluiD_Numerics_MeshFEMMove_StdUnSetup_hh

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
class StdUnSetup : public FEMMoveCom {
public:

  /**
   * Constructor.
   */
  explicit StdUnSetup(std::string name) : FEMMoveCom(name)
  {
  }

  /**
   * Destructor.
   */
  ~StdUnSetup()
  {
  }

  /**
   * Execute Processing actions
   */
  void execute();

}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshFEMMove

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshFEMMove_StdUnSetup_hh

