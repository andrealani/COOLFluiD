#ifndef COOLFluiD_Numerics_AeroCoef_NavierStokesSkinFrictionHeatFluxFS_hh
#define COOLFluiD_Numerics_AeroCoef_NavierStokesSkinFrictionHeatFluxFS_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataProcessingData.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace NavierStokes {
      class NavierStokesVarSet;
    }
  }

  namespace Numerics {
    
    namespace FluctSplit {
      class FluctuationSplitData;
      class InwardNormalsData;
    }
  }
  
  namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the skin friction and the heat flux for NavierStokes
 * simulations with
 * @see FluctuationSplit
 *
 * @author Andrea Lani
 * @author Thomas Wuilbaut
 *
 */
template <class UPDATEVAR>
class NavierStokesSkinFrictionHeatFluxFS : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  NavierStokesSkinFrictionHeatFluxFS(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~NavierStokesSkinFrictionHeatFluxFS();
  
  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

protected:

  /**
   * Execute on a set of dofs
   */
  void executeOnTrs();
  
  /**
   * Compute the required values
   */
  virtual void computeValues();

  /**
   * Compute the required values
   */
  virtual void computeExtraValues(Framework::GeometricEntity* currFace)
  {
  }
  
  /**
   * Compute the wall temperature and pressure
   */
  virtual void getWallTemperaturePressure(CFreal& TWall, CFreal& pWall);
  
  /**
   * Prepare Output file with the aerodynamic coef
   */
  virtual void prepareOutputFileAero();
  
  /**
   * Update the Output file with the aerodynamic coef
   */
  virtual void updateOutputFileAero();
 
  /**
   * Update the Output file with the wall values
   */
  virtual void prepareOutputFile(Framework::GeometricEntity* currFace);
    
  /**
   * Compute and output to screen the residuals of surface quantities of interest
   */
  virtual void computeSurfaceResiduals();
  
  /**
   * Write the surface files
   */
  virtual void writeSurfaceFiles();
  
  /**
   * Compute the value of tau at the wall
   */
  void computeTauWall(Framework::GeometricEntity* currFace);

  /// Update values to be printed and the corresponding residual
  void updateValuesMatAndResidual(CFuint iVar, CFuint index, CFreal value)
  {
    _valuesMatRes(iVar, index) = value - _valuesMat(iVar, index);
    _valuesMat(iVar, index) = value;
  }
  
protected:

  // the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  // the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;
  
  /// storage of the face normals
  Framework::DataSocketSink<FluctSplit::InwardNormalsData*> socket_normals;
  
  /// handle to the neighbor cell
  Framework::DataSocketSink<
    Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
  socket_faceNeighCell;
  
  // pointer to the data of the cell centered FVM method
  Common::SafePtr<FluctSplit::FluctuationSplitData> _fsData;
  
  /// diffusive variable set
  Common::SafePtr<Physics::NavierStokes::NavierStokesVarSet> _diffVar;
  
  /// update variable set
  Common::SafePtr<UPDATEVAR> _updateVar;
  
  // mapping between faceIDs and global index
  Common::CFMap<CFuint, CFuint> _mapTrsFaceToID;
  
  /// array storing the L2 norms of the values to write
  RealVector _valuesMatL2;

  /// array storing the L2 norms of the values to write
  RealVector _l2Norm;
  
  /// 2D array storing all values to write to file
  RealMatrix _valuesMat;
  
  /// 2D array storing the residuals of all values to write
  RealMatrix _valuesMatRes;
  
  // array of states involved in the flux computation
  std::vector<RealVector*> _states;
   
  // array of values (p, u, v, T, ...)
  RealMatrix _values;
  
  /// average cell values
  RealVector _avValues;
 
  /// average boundary face values
  Framework::State* _avBFaceValues;
 
  // temporary unit normal
  RealVector _unitNormal;

  // temporary unit normal
  RealVector _unitTangent;

  
  // temporary vibrational temperature
  RealVector _tempVib;
  
  // temporary coordinates of the cell center
  RealVector _coord;
  
  // physical data array
  RealVector _pData;
    
  //temporary value of density
  CFreal _rhoWall;

  //temporary value of dyn viscosity
  CFreal _muWall;
  
  //tau at the wall in the 2D case
  CFreal _tau;
  
  //skinfriction in the 2D case
  CFreal _skinFriction;
  
  //tau at the wall in the 2D case
  RealVector _frictionDrag;

  //Lift, Drag
  CFreal _lift;
  CFreal _drag;
  
  /// append time to the files names
  bool _appendTime;
  
  /// append the iteration number to the files names
  bool _appendIter;
    
  // Storage for choosing when to save the wall values file
  CFuint _saveRate;
  
  // name of the surface convergence file
  std::string _outputFileConv;
  
  // name of the output file
  std::string _outputFile;

  // name of the output file
  std::string _outputFileAero;

  //freestream values
  CFreal _uInf;
  CFreal _pInf;
  CFreal _TInf;

  // ID of temperature in gradient vars
  CFuint _TID;

  // IDs of velocity component Vx in gradient vars
  CFuint _UID;

  // IDs of velocity component Vy in gradient vars
  CFuint _VID;

  // IDs of velocity component Vz in gradient vars
  CFuint _WID;

}; // end of class NavierStokesSkinFrictionHeatFluxFS

//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokesSkinFrictionHeatFluxFS.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_AeroCoef_NavierStokesSkinFrictionHeatFluxFS_hh
