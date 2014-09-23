#ifndef COOLFluiD_Numerics_SubSystemCoupler_StdPostProcess_hh
#define COOLFluiD_Numerics_SubSystemCoupler_StdPostProcess_hh

//////////////////////////////////////////////////////////////////////////////

#include "SubSysCouplerData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to Domain to be executed in order to setup the MeshData.
   *
   * @author Thomas Wuilbaut
   *
   */
class StdPostProcess : public CouplerCom {
public:

  /**
   * Constructor.
   */
  explicit StdPostProcess(std::string name) : CouplerCom(name)
  {
  }

  /**
   * Destructor.
   */
  ~StdPostProcess()
  {
  }

protected:

  /**
   * Execute Processing actions
   */
  void executeOnTrs();

}; // class StdPostProcess

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_StdPostProcess_hh

