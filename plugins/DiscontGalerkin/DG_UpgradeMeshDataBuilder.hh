#ifndef COOLFluiD_DiscontGalerkin_DG_UpgradeMeshDataBuilder_hh
#define COOLFluiD_DiscontGalerkin_DG_UpgradeMeshDataBuilder_hh

//////////////////////////////////////////////////////////////////////////////

#include "DiscontGalerkin/DG_MeshDataBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

 namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

/**
* This class builds data inside MeshData.
* It assumes a Finite Volume mesh has been read and it upgrades
* the elements to DiscontGalerkin elements
*/
class DG_UpgradeMeshDataBuilder : public DiscontGalerkin::DG_MeshDataBuilder {

public:// static functions

 /**
 * Defines the Config Option's of this class
 * @param options a OptionList where to add the Option's
 */
 static void defineConfigOptions(Config::OptionList& options);

public:// functions

 /**
 * Constructor
 *
 * @param name name of the builder used for configuration
 */
 DG_UpgradeMeshDataBuilder(const std::string& name);

 /**
 * Destructor
 */
 ~DG_UpgradeMeshDataBuilder();

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
 CFuint getNbrOfStatesInType(CFGeoShape::Type geoShape, CFPolyOrder::Type polyOrder);

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

protected:// functions

 /**
 * Recreate the cell-states connectivity
 * for the new DiscontGalerkin elements
 */
 virtual void upgradeStateConnectivity();

 /**
 * Recreates the states to match the new DiscontGalerkin elements
 */
 virtual void recreateStates();

protected:// data

 /// spectral finite volume polynomial order
 CFPolyOrder::Type m_PolyOrder;

 /// string holding the spectral finite volume polynomial order
 std::string m_PolyOrderStr;

}; // end of class DG_UpgradeMeshDataBuilder

//////////////////////////////////////////////////////////////////////////////

 } // namespace DiscontGalerkin
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_DiscontGalerkin_DG_UpgradeMeshDataBuilder_hh

