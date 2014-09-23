#ifndef COOLFluiD_Numerics_MeshAdapterSpringAnalogy_ImposedValvePrepare_hh
#define COOLFluiD_Numerics_MeshAdapterSpringAnalogy_ImposedValvePrepare_hh

//////////////////////////////////////////////////////////////////////////////

#include "BasePrepare.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshAdapterSpringAnalogy {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a NumericalCommand action to be
   * sent to be executed in order to move the boundaries of the mesh in the
   * case of a moving valve
   *
   * @author Thomas Wuilbaut
   *
   */

class ImposedValvePrepare : public BasePrepare {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  ImposedValvePrepare(const std::string& name);

  /**
   * Destructor.
   */
  ~ImposedValvePrepare()
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

  /// name of the trs of the valve
  std::string _trsNameValve;

  /// name of the trs of the symmetry plane in cylinder
  std::string _trsNameSymmetry;

  /// Valve translation distance
  CFreal _translation;

  /// Domain maximum X coordinate
  CFreal _maxXCoord;

  /// Domain minimum X coordinate
  CFreal _minXCoord;

  /// Valve minimum X coordinate where Y is constant
  CFreal _minXCoord_Yconst;

}; // class ImposedValvePrepare

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshAdapterSpringAnalogy

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshAdapterSpringAnalogy_ImposedValvePrepare_hh

