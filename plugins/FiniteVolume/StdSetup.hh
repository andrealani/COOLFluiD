#ifndef COOLFluiD_Numerics_FiniteVolume_StdSetup_hh
#define COOLFluiD_Numerics_FiniteVolume_StdSetup_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseSetupFVMCC.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a command to be executed during the set up
 * of a standard Finite Volume Method.
 *
 * @author Andrea Lani
 */
class StdSetup : public Framework::BaseSetupFVMCC<CellCenterFVMCom> {
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  explicit StdSetup(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~StdSetup();
  
  /**
   * Configure the command
   */
  virtual void configure ( Config::ConfigArgs& args );
  
  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Execute Processing actions
   */
  virtual void execute();

protected: // helper method
  
  /**
   * Compute all 3D geometric data
   */
  void compute3DGeometricData();
  
  /**
   * Check if the centroid is inside (only working for hexahedra for now)
   */
  bool isCentroidInside(CFuint cellID, 
			const std::vector<Framework::GeometricEntity*>& facesInCell, 
			const RealVector& xc); 
  
  /**
   * Check if the centroid is inside (only working for hexahedra for now)
   */
  void isCentroidInside();
  
  /**
   * Assign the face nodes
   */
  void assignFaceNodes(Framework::GeometricEntity *const cell);
  
  /**
   * Compute the face geometry
   */
  void computeFaceGeometry(const Framework::GeometricEntity *const face, 
			   RealVector& normalPtr, 
			   RealVector& xcFacePtr);
  /**
   * Compute the face geometry
   */
  void computeFaceGeometryIter(const Framework::GeometricEntity *const face, 
			       RealVector& normalPtr, 
			       RealVector& xcFacePtr);
  
  /**
   * Compute the cell geometry
   */
  void computeCellGeometry3D(CFuint iElem, 
			     const Framework::GeometricEntity *const cell,
			     RealVector& xcFacePtr,
			     CFuint& countNegativeVolumes,
			     CFuint& count);
  
  /**
   * Compute the cell geometry
   */
  void computeCellGeometry3DIter(CFuint iElem, 
				 const Framework::GeometricEntity *const cell, 
				 RealVector& xcFacePtr,
				 CFuint& countNegativeVolumes, 
				 CFuint& count);
  
protected:
    
  // flags telling if diffusion terms have to be computed for a given cell
  Framework::DataSocketSource<CFreal> socket_activeDiffusion;
  
  // type of the stencil
  std::string _stencilType;
  
  /// initial limiter socket name
  std::string _limiterSocketName;
  
  // temporary array
  RealVector _v1;
  
  // temporary array
  RealVector _v2; 
  
  // temporary array
  RealVector _v3;

  // temporary array
  RealVector _xcCellApprox;
      
}; // class Setup

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_StdSetup_hh

