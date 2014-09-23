#ifndef COOLFluiD_Numerics_SubSystemCoupler_FVMCCPreProcessWrite_hh
#define COOLFluiD_Numerics_SubSystemCoupler_FVMCCPreProcessWrite_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/FaceTrsGeoBuilder.hh"
#include "SubSystemCoupler/StdPreProcessWrite.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a NumericalCommand action to be
 * to be executed in order to setup the SubSystemCoupler Method
 *
 * @author Thomas Wuilbaut
 *
 */
class FVMCCPreProcessWrite : public StdPreProcessWrite {
public:

  /**
   * Constructor.
   */
  explicit FVMCCPreProcessWrite(const std::string& name);

  /**
   * Destructor.
   */
  ~FVMCCPreProcessWrite();

protected:
  /**
   * Execute Processing actions
   */
  void executeOnTrs();

  /**
   * Determine and fill the datahandle of coordinates
   * in case the transfer is done at the nodes
   * @param socketName name of the socket to store the coordinates of the nodes
   */
  void fillGhostCoordDataHandle(const std::string& socketName);

  /**
   * Determine and fill the datahandle of coordinates
   * in case the transfer is done at the nodes
   * @param socketName name of the socket to store the coordinates of the nodes
   */
  void fillNodalCoordDataHandle(const std::string& socketName);

  /**
   * Determine and fill the datahandle of coordinates
   * in case the transfer is done at the nodes
   * @param socketName name of the socket to store the coordinates of the nodes
   */
  void fillStatesCoordDataHandle(const std::string& socketName);

  /**
   * In case the type of values to be passed is Nodal,
   * Create a datahandle with a vector of pointers to
   * the Geometric entities containing that node
   * @param socketName name of the socket that contains the NodeToFace connectivity
   */
  virtual void setNodeToFaceConnectivity(const std::string& socketName);

protected: // data

}; // class FVMCCPreProcessWrite

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_FVMCCPreProcessWrite_hh

