#ifndef COOLFluiD_Numerics_SubSystemCoupler_FVMCCNewtonMeshMatcherWrite_hh
#define COOLFluiD_Numerics_SubSystemCoupler_FVMCCNewtonMeshMatcherWrite_hh

//////////////////////////////////////////////////////////////////////////////

#include "SubSystemCoupler/NewtonMeshMatcherWrite.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * sent to Domain to be executed in order to set the
 * match between meshes
 *
 * @author Thomas Wuilbaut
 */
class FVMCCNewtonMeshMatcherWrite : public NewtonMeshMatcherWrite {
public:

  /**
   * Constructor.
   */
  explicit FVMCCNewtonMeshMatcherWrite(const std::string& name);

  /**
   * Destructor.
   */
  ~FVMCCNewtonMeshMatcherWrite();

protected:

  /**
   * Executes the command.
   */
  virtual void executeWrite(const CFuint iProc);


private: // functions

  /**
   * Preselect the faces on which to apply the projection algorithm
   */
  virtual void facePreSelection();

  /**
   * Pairs the projected point with the closest face
   */
  virtual void nodeToElementPairing(SubSysCouplerData::GeoEntityIdx& matchingFace);

}; // class FVMCCNewtonMeshMatcherWrite

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_FVMCCNewtonMeshMatcherWrite_hh

