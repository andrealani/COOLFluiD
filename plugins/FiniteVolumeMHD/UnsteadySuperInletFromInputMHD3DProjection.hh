#ifndef COOLFluiD_Numerics_FiniteVolume_UnsteadySuperInletFromInputMHD3DProjection_hh
#define COOLFluiD_Numerics_FiniteVolume_UnsteadySuperInletFromInputMHD3DProjection_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "FiniteVolume/SuperInlet.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {
    namespace MHD {
      class MHD3DProjectionVarSet;
    }
  }

  namespace Numerics {
    
    namespace FiniteVolume {
      
      namespace FiniteVolumeMHD {
	
//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents an unsteady superfast inlet command based on an input file
   *
   * @author Thomas Wuilbaut
   * @author Mehmet Sarp Yalim
   *
   */
class UnsteadySuperInletFromInputMHD3DProjection : public SuperInlet {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  UnsteadySuperInletFromInputMHD3DProjection(const std::string& name);

  /**
   * Default destructor
   */
  ~UnsteadySuperInletFromInputMHD3DProjection();

  /** 
   * Set the preProcesses connectivity between faces belonging to different process 
   * 
   */ 
  void preProcess(); 

  /**
   * Read the inlet solar wind data from ACE satellite from the input file
   */
  void readInputFile();

  /** 
   * Read the inlet solar wind data from ACE satellite from the input file 
   */ 
  void readInputFileNoInterpolation(); 
  
  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

  /**
   * Configures the command.
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Apply boundary condition on the given face
   */
  void setGhostState(Framework::GeometricEntity *const face);

protected: // data and methods

  /**
   * Construct the file path from which to read the solar wind data
   */
  boost::filesystem::path constructFilename();

private:
  
  /// corresponding variable set
  Common::SafePtr<Physics::MHD::MHD3DProjectionVarSet> _varSet;
  
  /// storage for the time from the input file
  RealVector _t;
  
  /// storage for the proton density from the input file
  RealVector _rho;
  
  /// storage for the x-velocity of the solar wind from the input file
  RealVector _u;

  /// storage for the y-velocity of the solar wind from the input file
  RealVector _v;

  /// storage for the z-velocity of the solar wind from the input file
  RealVector _w;

  /// storage for the x-component of the magnetic field of the solar wind from the input file
  RealVector _bx;

  /// storage for the y-component of the magnetic field of the solar wind from the input file
  RealVector _by;

  /// storage for the z-component of the magnetic field of the solar wind from the input file
  RealVector _bz;

  /// physical model data
  RealVector _dataInnerState;

  /// storage for the pressure of the solar wind from the input file
  RealVector _p;

  /// Name of input file from where to read the ACE solar wind data
  std::string _nameInputFile;
  
  /// Begin time of the overall unsteady simulation that uses ACE data input file in case of restart 
  CFreal _beginTime;
  
  /// Maximum time of the unsteady simulation that uses ACE data input file 
  CFreal _maxTime;
  
  /// flag to make sure that the input file for the solar wind data from ACE satellite is read only once
  CFint _flagReadInputFile;

  /// flag telling to use time interpolation while reading the input file
  bool _useTimeInterpolation;

}; // end of class UnsteadySuperInletFromInputMHD3DProjection

//////////////////////////////////////////////////////////////////////////////

 } // namespace FiniteVolumeMHD

  } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_UnsteadySuperInletFromInputMHD3DProjection_hh
