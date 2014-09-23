#ifndef COOLFluiD_Physics_StructMech_Materials_hh
#define COOLFluiD_Physics_StructMech_Materials_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/CFMap.hh"
#include "MathTools/RealMatrix.hh"
#include "MathTools/MatrixInverterT.hh"
#include "StructMech/MaterialData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a library for Materials characteristics.
 *
 * @author Thomas Wuilbaut
 *
 */
class Materials{
public:

  /**
   * Constructor without arguments
   */
  Materials();

  /**
   * Default destructor
   */
  ~Materials();

  /**
   * Get the Physical data
   */
  MaterialData getMaterialData(const std::string name)
  {
    cf_assert(_data.exists(name));

    return _data.find(name);
  }

private:
  /**
   * Define the Material Properties
   */
  void addSteel();
  void addAluminium();
  void addTitanium();
  void addComposites();

protected:

  /// Material Data Map
  Common::CFMap<std::string,MaterialData> _data;

  /// Temporary Material Data
  MaterialData _temp;
  /// matrix inverter size 3
  MathTools::MatrixInverterT<3> m_inverter3;

}; // end of class Materials

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_StructMech_Materials_hh
