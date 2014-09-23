#ifndef COOLFluiD_SpectralFV_MeshUpgradeBuilder_hh
#define COOLFluiD_SpectralFV_MeshUpgradeBuilder_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFV/SpectralFVBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class builds data inside MeshData.
 * It assumes a Finite Volume mesh has been read and it upgrades
 * the elements to SpectralFV elements
 */
class MeshUpgradeBuilder : public SpectralFV::SpectralFVBuilder {

public: // static functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

public: // functions

  /**
   * Constructor
   *
   * @param name name of the builder used for configuration
   */
  MeshUpgradeBuilder(const std::string& name);

  /**
   * Destructor
   */
  ~MeshUpgradeBuilder();

  /**
   * Releases temporary memory used in building the mesh
   */
  virtual void releaseMemory();

  /// Configures this object from the supplied arguments
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the number of control volumes in a spectral volume of given shape and order
   * @param geoShape the shape of the spectral volume
   * @param polyOrder the order of the polynomial interpolation in the spectral volume
   * @return the number of control volumes in the given spectral volume type
   */
  CFuint getNbrOfStatesInSVType(CFGeoShape::Type geoShape, CFPolyOrder::Type polyOrder);

  /**
   * Creates all the TopologicalRegionSet's in the MeshData
   * This function overides the one in MeshDataBuilder
   */
  virtual void createTopologicalRegionSets();

  /**
   * Computes the data for the element types
   * This function overides the one in MeshDataBuilder
   */
  virtual void computeGeoTypeInfo();

protected: // functions

  /**
   * Recreate the cell-states connectivity
   * for the new SpectralFV elements
   */
  virtual void upgradeStateConnectivity();

  /**
   * Recreates the states to match the new SpectralFV elements
   */
  virtual void recreateStates();

protected: // data

  /// spectral finite volume polynomial order
  CFPolyOrder::Type m_svPolyOrder;

  /// string holding the spectral finite volume polynomial order
  std::string m_svPolyOrderStr;

};  // end of class MeshUpgradeBuilder

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_SpectralFV_MeshUpgradeBuilder_hh

