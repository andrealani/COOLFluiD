#ifndef COOLFluiD_Numerics_MeshAdapterSpringAnalogy_ImposedMovementPrepare_hh
#define COOLFluiD_Numerics_MeshAdapterSpringAnalogy_ImposedMovementPrepare_hh

//////////////////////////////////////////////////////////////////////////////

#include "BasePrepare.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshAdapterSpringAnalogy {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to Domain to be executed in order to setup the MeshData.
   *
   * @author Thomas Wuilbaut
   *
   */

class ImposedMovementPrepare : public BasePrepare {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  ImposedMovementPrepare(const std::string& name);

  /**
   * Destructor.
   */
  ~ImposedMovementPrepare()
  {
  }

  /**
   * Configures this object with supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

private: // functions

  /**
   * Move the boundary nodes
   */
  virtual void moveBoundaries();

private: // data

  /// name of the trs
  std::string _trsName;

  /// rotation angle
  CFreal _rotationAngle;

  /// rotation center for config
  std::vector<CFreal> _rotationCenterConf;

  /// rotation center
  RealVector _rotationCenter;

}; // class ImposedMovementPrepare

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshAdapterSpringAnalogy

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshAdapterSpringAnalogy_ImposedMovementPrepare_hh

