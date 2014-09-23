#ifndef COOLFluiD_Numerics_AeroCoef_NavierStokesConsComputeAeroSpectralFDCompact_hh
#define COOLFluiD_Numerics_AeroCoef_NavierStokesConsComputeAeroSpectralFDCompact_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "MathTools/FunctionParser.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/DynamicDataSocketSet.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/FaceTrsGeoBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

    namespace Framework
    {
      class FaceToCellGEBuilder;
    }
    namespace Physics
    {
      namespace NavierStokes
      {
        class EulerVarSet;
        class NavierStokesVarSet;
      }
    }
    namespace SpectralFD
    {
      class CompactBndFaceTermComputer;
      class BCStateComputer;
    }

//////////////////////////////////////////////////////////////////////////////

    namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the Wall values and aerodynamic coefficients
 * for NavierStokes2DCons with spectral difference
 *
 * @author Kris Van den Abeele
 *
 */
class NavierStokesConsComputeAeroSpectralFDCompact : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  NavierStokesConsComputeAeroSpectralFDCompact(const std::string& name);

  /**
   * Default destructor
   */
  ~NavierStokesConsComputeAeroSpectralFDCompact();

  /**
   * Configure the command
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command
   */
  void unsetup();

protected:

  /**
   * Execute on a set of dofs
   */
  void executeOnTrs();

  /**
   * Compute the values at the wall and write them to file
   */
  void computeWall();

  /**
   * Compute the aerodynamic coefficients and write them to file
   */
  void computeAero();

  /**
   * Open the Output File and Write the header
   */
  void prepareOutputFileWall();

  /**
   * Open the Output File and Write the header
   */
  void prepareOutputFileAero();

  /**
   * Write the aerodynamic coefficients to file
   */
  void updateOutputFileAero();

  /// compute the required data in the face term for the current face
  void computeFaceData();

private:

  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> m_sockets;

  /// builder of faces
  Common::SafePtr<Framework::GeometricEntityPool<Framework::FaceToCellGEBuilder> > m_faceBuilder;

  /// SpectralFD strategy that computes the boundary face terms
  Common::SafePtr< SpectralFD::CompactBndFaceTermComputer > m_bndFaceTermComputer;

  /// the BCStateComputer for this BC
  Common::SafePtr< SpectralFD::BCStateComputer > m_bcStateComputer;

  /// Update variable set
  Common::SafePtr<Physics::NavierStokes::EulerVarSet> m_updateVarSet;

  /// Diffusive variable set
  Common::SafePtr<Physics::NavierStokes::NavierStokesVarSet> m_diffusiveVarSet;

  /// variable for current face
  Framework::GeometricEntity* m_face;

  /// variable for current internal cell
  Framework::GeometricEntity* m_intCell;

  /// variable for current face orientation
  CFuint m_orient;

  /// the states in the neighbouring cell
  std::vector< Framework::State* >* m_cellStates;

  /// Storage for choosing when to save the wall values file
  CFuint m_saveRateWall;

  /// Storage for choosing when to save the aero coef file
  CFuint m_saveRateAero;

  /// Output File for Wall Values
  Common::SelfRegistPtr<Environment::FileHandlerOutput> m_fileWall;

  /// Name of Output File where to write the coeficients.
  std::string m_nameOutputFileWall;

  /// Name of Output File where to write the coeficients.
  std::string m_nameOutputFileAero;

  //// Variables to store the wall values + lift/drag
  /// Storage for lift coeficient on trs
  CFreal m_lift;

  /// Storage for drag coeficient on trs
  CFreal m_drag;

  /// Storage for friction drag coeficient on trs
  CFreal m_fricDrag;

  /// force coefficient
  RealVector m_forceCoef;

  /// friction force coefficient
  RealVector m_fricForceCoef;

  /// velocity at infinity
  CFreal m_uInf;

  /// density at infinity
  CFreal m_rhoInf;

  /// pressure at infinity
  CFreal m_pInf;

  /// flow direction
  std::vector< CFreal > m_flowDir;

  /// projected surface
  CFreal m_projSurf;

  ///flag for appending iteration
  bool m_appendIter;

  ///flag for appending time
  bool m_appendTime;

  /// number of equations in the physical model
  CFuint m_nbrEqs;

  /// variable containing current gradients
  std::vector< RealVector*  > m_grads;

}; /// end of class NavierStokesConsComputeAeroSpectralFDCompact

//////////////////////////////////////////////////////////////////////////////

    } /// namespace AeroCoef

} /// namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif /// COOLFluiD_Numerics_AeroCoef_NavierStokesConsComputeAeroSpectralFDCompact_hh
