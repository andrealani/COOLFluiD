#ifndef COOLFluiD_Numerics_MeshFEMMove_Valve3DCyclePrepare_hh
#define COOLFluiD_Numerics_MeshFEMMove_Valve3DCyclePrepare_hh

//////////////////////////////////////////////////////////////////////////////

#include "BasePrepare.hh"
#include "Framework/DynamicDataSocketSet.hh"

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

class Valve3DCyclePrepare : public BasePrepare {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  Valve3DCyclePrepare(const std::string& name);

  /**
   * Destructor.
   */
  ~Valve3DCyclePrepare()
  {
  }

  /**
   * Configures this object with supplied arguments.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set up the data
   */
  void setup();

  /**
   * Returns the DataSocket's that this command needs as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

private: // functions

  /**
   * Move the boundary nodes
   */
  virtual void moveBoundaries();

private: // data

  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> _sockets;

  /// name of the trs of the valve
  std::string _trsNameValve;

  /// Valve translation distance
  CFreal _translation;

  /// Maximum Valve opening
  CFreal _maxOpening;
  CFreal _minOpening;

  /// Engine Rotation Speed
  CFreal _rpm;

  /// Domain maximum X coordinate
  CFreal _maxXCoord;

  /// Domain minimum X coordinate
  CFreal _minXCoord;

  /// Valve minimum X coordinate where Y is constant
  CFreal _minXCoord_Yconst;

}; // class Valve3DCyclePrepare

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshFEMMove

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_MeshFEMMove_Valve3DCyclePrepare_hh

