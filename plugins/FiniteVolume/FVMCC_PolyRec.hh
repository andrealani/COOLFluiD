#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_PolyRec_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_PolyRec_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/PolyReconstructor.hh"
#include "Framework/GeometricEntityPool.hh"
#include "Framework/BaseDataSocketSink.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/VectorialFunction.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {
      class CellCenterFVMData;
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class implements a constant polynomial reconstructor for FVM
 *
 * @author Andrea Lani
 *
 */
class FVMCC_PolyRec : public Framework::PolyReconstructor<CellCenterFVMData> {
public:
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
  
  /**
   * Constructor
   */
  FVMCC_PolyRec(const std::string& name);
  
  /**
   * Default destructor
   */
  virtual ~FVMCC_PolyRec();
  
  /**
   * Set private data that will be used during the computation
   */
  virtual void setup();
  
  /**
   * Unsetup private data that will be used during the computation
   */
  virtual void unsetup();

  /**
   * Update the weights when nodes are moving
   */
  virtual void updateWeights();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
 
  /**
   * Compute the gradients
   */
  virtual void computeGradients() = 0;
  
  /// Get the current left state
  Framework::State& getCurrLeftState()
  {
    return *_extrapValues[0];
  }
  
  /// Get the current right state
  Framework::State& getCurrRightState()
  {
    return *_extrapValues[1];
  }
  
  /**
   * Return the nb of quadrature points per face
   */
  CFuint nbQPoints() const {return 1;}
  
  /// Get the Geometric Shape Functions at current quadrature point
  const RealVector& getCurrentGeoShapeFunction(Framework::GeometricEntity *const face);
  
  /**
   * Prepare the reconstruction step
   */
  virtual void prepareReconstruction();
  
protected: // helper functions
  
  /**
   * Compute the limiters for a given face
   */
  void computeFaceLimiter(Framework::GeometricEntity *const face);
  
  /**
   * Compute the mid point
   */
  void computeMidPoint(const std::vector<Framework::Node*>& faceNodes, 
		       Framework::Node& faceMidCoord) 
  {
    faceMidCoord = 0.0;
    const CFuint nbFaceNodes = faceNodes.size();
    const CFreal ovNbFaceNodes = 1./nbFaceNodes;
    for (CFuint i = 0; i < nbFaceNodes; ++i) {
      faceMidCoord += (*faceNodes[i])*ovNbFaceNodes;
    }
  }  
  
  /// Allocate reconstruction data needed for the flux evaluation
  virtual void allocateReconstructionData();
  
  /**
   * Extrapolate the solution in the face quadrature points
   */
  void baseExtrapolateImpl(Framework::GeometricEntity* const face);
  
  /**
   * Constantly extrapolate the solution in the face quadrature points
   * This ensures a default exrapolation
   */
  void constantRecImpl(Framework::GeometricEntity* const face)
  {
    getValues(LEFT).copyData(*face->getState(LEFT));
    getValues(RIGHT).copyData(*face->getState(RIGHT));
    getBackupValues(LEFT) = getValues(LEFT);
    getBackupValues(RIGHT) = getValues(RIGHT);
  }
    
  /**
   * Constantly extrapolate the solution in the face quadrature points
   * This ensures a default exrapolation
   */
  void constantRecImpl(Framework::GeometricEntity* const face,
		       CFuint iVar, CFuint leftOrRight)
  {
    // for each quadrature point, set
    // reconstructed variable = same variable in the current state
    copyBackupValues();
    // set the new values
    getValues(leftOrRight)[iVar]   = (*face->getState(leftOrRight))[iVar];
  }
  
protected:
  
  /// flags for cells
  Framework::DataSocketSink<bool> socket_cellFlag;

  /// socket for nodes
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// socket for states
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// socket for ghost states
  Framework::DataSocketSink< Framework::State*> socket_gstates;
   
  /// socket for face normals
  Framework::DataSocketSink< CFreal> socket_normals;
  
  /// limiter is null
  bool _isLimiterNull;
  
  /// storage for the temporary extrapolated coordinates
  /// for all the faces of the current cell
  std::vector<std::vector<Framework::Node*> > _quadPointCoord;
  
  /// temporary limiter value
  RealVector _tmpLimiter; 
    
  /// gradient coefficient
  RealVector _gradientCoeff;
    
  /// the VectorialFunction to use
  Framework::VectorialFunction _vFunction;
  
  /// a vector of string to hold the functions
  std::vector<std::string> _functions;
  
  /// a vector of string to hold the functions
  std::vector<std::string> _vars;

  /// flag to stop the limiting
  CFuint _stopLimiting;
  
}; // end of class FVMCC_PolyRec

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_PolyRec_hh
