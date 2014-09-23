#ifndef COOLFluiD_Numerics_SubSystemCoupler_StructPreProcessWrite_hh
#define COOLFluiD_Numerics_SubSystemCoupler_StructPreProcessWrite_hh

//////////////////////////////////////////////////////////////////////////////

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
class StructPreProcessWrite : public StdPreProcessWrite {
public:

  /**
   * Constructor.
   */
  explicit StructPreProcessWrite(const std::string& name);

  /**
   * Destructor.
   */
  ~StructPreProcessWrite();

protected:

  /**
   * Determine and fill the datahandle of coordinates
   * in case the transfer is done at the nodes
   * @param socketName name of the socket to store the coordinates of the nodes
   */
  void fillNodalCoordDataHandle(const std::string& socketName);


}; // class StructPreProcessWrite

//////////////////////////////////////////////////////////////////////////////

    } // namespace SubSystemCoupler

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_SubSystemCoupler_StructPreProcessWrite_hh

