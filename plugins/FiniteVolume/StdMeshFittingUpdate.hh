#ifndef COOLFluiD_Numerics_FiniteVolume_StdMeshFittingUpdate_hh
#define COOLFluiD_Numerics_FiniteVolume_StdMeshFittingUpdate_hh

//////////////////////////////////////////////////////////////////////////////

#include "CellCenterFVMData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a command to be executed after
   * the mesh has been updated
   */

//////////////////////////////////////////////////////////////////////////////

class StdMeshFittingUpdate : public CellCenterFVMCom {
public:

  /**
   * Constructor.
   */
  explicit StdMeshFittingUpdate(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~StdMeshFittingUpdate();
  
  /**
   * Execute Processing actions
   */
  virtual void execute();
  
}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_StdMeshFittingUpdate_hh

