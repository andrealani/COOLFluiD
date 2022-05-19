#ifndef COOLFluiD_Numerics_FluxReconstructionMethod_BCSuperInletProjMHD_hh
#define COOLFluiD_Numerics_FluxReconstructionMethod_BCSuperInletProjMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/BaseMethodStrategyProvider.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "Framework/FaceTrsGeoBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
    
  namespace Framework
  {
    class FaceToCellGEBuilder;
  }
    
  namespace Physics {
    namespace MHD {
      class MHD3DProjectionVarSet;
    }
  }

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents an inlet boundary condition
 * for the 3D MHD equations.
 *
 * @author Ray Vandenhoeck
 */
class BCSuperInletProjMHD : public BCStateComputer {

public:  // methods

  /// Constructor
  BCSuperInletProjMHD(const std::string& name);

  /// Destructor
  ~BCSuperInletProjMHD();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "BCSuperInletProjMHD";
  }
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /// Set up private data and data
  void setup();

  /**
   * Sets the ghost states in all the boundary points (depends on the boundary condition type)
   */
  void computeGhostStates(const std::vector< Framework::State* >& intStates,
                          std::vector< Framework::State* >& ghostStates,
                          const std::vector< RealVector >& normals,
                          const std::vector< RealVector >& coords);

  /**
   * Sets the ghost gradients in all the boundary points (depends on the boundary condition type)
   */
   void computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                              std::vector< std::vector< RealVector* > >& ghostGrads,
                              const std::vector< RealVector >& normals,
                              const std::vector< RealVector >& coords);
   
   /**
   * Set the B boundary values
   */
  virtual void preProcess();


protected: // data
  
  /// map TRS name -> initial solution array that will be used as BC value
  Common::CFMap<std::string, RealVector*> m_initialSolutionMap;
  
  ///bnd density value
  CFreal m_rhoBC;
  
  ///bnd p value
  CFreal m_pBC;
  
  /// array specifying IDs of initial solution components that will be used as BC value
  std::vector<CFuint> m_initialSolutionIDs;
  
  /// builder of faces
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> > m_faceBuilder;
  
  /// variable for current internal cell
  Framework::GeometricEntity* m_intCell;
  
  /// the states in the neighbouring cell
  std::vector< Framework::State* >* m_cellStates;
  
  /// number of flux pnts on a face
  CFuint m_nbrFaceFlxPnts;
  
  /// current face
  Framework::GeometricEntity* m_currFace;
  
  /// variable for current face orientation
  CFuint m_orient;
  
  /// extrapolated states in the flux points of the cell
  std::vector< Framework::State* > m_cellStatesFlxPnt;
  
  /// coefs to extrapolate the states to the flx pnts
  Common::SafePtr< std::vector< std::vector< CFreal > > > m_solPolyValsAtFlxPnts;
  
  /// flx pnt - face connectivity
  Common::SafePtr< std::vector< std::vector< CFuint > > > m_faceFlxPntConn;
  
  /// TRS in which this BC is defined
  Common::SafePtr< Framework::TopologicalRegionSet > m_thisTRS;

}; // class BCSuperInletProjMHD

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_FluxReconstructionMethod_BCSuperInletProjMHD_hh

