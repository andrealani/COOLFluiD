#ifndef COOLFluiD_Numerics_AeroCoef_NavierStokesBLExtractionCC_hh
#define COOLFluiD_Numerics_AeroCoef_NavierStokesBLExtractionCC_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesSkinFrictionHeatFluxCC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class extracts a BL from a 2D body
 *
 * @author Thomas Wuilbaut
 *
 */
class NavierStokesBLExtractionCC : public NavierStokesSkinFrictionHeatFluxCC {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  NavierStokesBLExtractionCC(const std::string& name);

  /**
   * Default destructor
   */
  ~NavierStokesBLExtractionCC();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

protected:

  /**
   * Execute on a set of dofs
   */
  void executeOnTrs();

  /**
   * Where the actual extraction along a normal to a face takes place
   */
  virtual void computeExtraValues();
  
  /**
   * Extract Boundary Layer Profile
   */
  void extractBoundaryLayerProfile(Framework::GeometricEntity* currFace);

  /**
   * Find Stagnation Point
   */
  void findStagnationPoint();

  /**
   * Extract boundary layer profiles along the TRS
   */
  void extractBLalongProfile();

  /**
   * Build Node Face connectivity
   */
  void buildNodeFaceConnectivity();

protected:

  //name of output file for the BL profile
  std::string _outputFileBL;

  //coordinates for the extraction of BL profile
  std::vector<CFreal> _extractCoord;
  RealVector _initExtract;

  //thickness of the BL
  CFreal _BLThickness;

  //Node index (in the trs) of the stagnation point and of the trailing edge
  CFuint _trailingEdgeIdx;
  CFuint _stagnationPointIdx;

  // node to face connectivity
  std::vector< std::vector<CFuint> > _nodeToFaceConnectivity;

  //Build the nodes LocalID to TrsIdx map.
  Common::CFMap<CFuint, CFuint> _nodesLocalID2TrsIdx;

  //distance between stagnation point and the next node
  CFreal _distanceToStagnation;

  //flag to tell if need to extract BL along a profile
  bool _extractBLalongProfile;

  //Flag to check if the profile to extract could be found
  bool _initCoordFound;

  //tolerance for finding the initial coord on a face
  CFreal _tolerance;

}; // end of class NavierStokesBLExtractionCC

//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_AeroCoef_NavierStokesBLExtractionCC_hh
