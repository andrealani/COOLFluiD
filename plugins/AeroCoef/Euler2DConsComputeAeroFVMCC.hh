#ifndef COOLFluiD_Numerics_AeroCoef_Euler2DConsComputeAeroFVMCC_hh
#define COOLFluiD_Numerics_AeroCoef_Euler2DConsComputeAeroFVMCC_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "MathTools/FunctionParser.hh"
#include "Framework/DataProcessingData.hh"
#include "Framework/DynamicDataSocketSet.hh"
#include "Framework/DataSocketSink.hh"
#include "Framework/FaceTrsGeoBuilder.hh"
#include "NavierStokes/Euler2DCons.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the Wall values and aerodynamic coefficients
 * for Euler2DCons
 *
 * @author Thomas Wuilbaut
 *
 */
class Euler2DConsComputeAeroFVMCC : public Framework::DataProcessingCom {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  Euler2DConsComputeAeroFVMCC(const std::string& name);

  /**
   * Default destructor
   */
  ~Euler2DConsComputeAeroFVMCC();

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
  void computeAeroFVMCC();
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

private:

  /**
   * Set the functions for Alpha angle
   * @throw Common::ParserException if the expression is senseless
   * @throw BadValueException if the expression defines more than one functions for alpha
   */
  void setFunction();

private:

  /// builder for Cell Centered FVM TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::FaceTrsGeoBuilder> m_faceTrsGeoBuilder;

  /// the dynamic sockets in this Command
  Framework::DynamicDataSocketSet<> m_sockets;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::Node* , Framework::GLOBAL > socket_nodes;

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// the socket to the data handle of the nodal state's
  Framework::DataSocketSink<RealVector> socket_nstates;

  /// the socket to the data handle of the nodal state's
  Framework::DataSocketSink<Framework::State*> socket_gstates;

  /// Update variable set
  Common::SafePtr<Physics::NavierStokes::EulerVarSet> m_updateVarSet;

  /// physical model data
  RealVector m_dataState;

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

  /// Alpha function of time
  MathTools::FunctionParser m_functionParser;

  /// a vector of string to hold the functions
  std::string m_function;

  /// a vector of string to hold the functions
  std::string m_vars;

  //// Variables to store the wall values + lift/drag
  /// pressure
  CFreal m_p;

  /// pressure coefficient
  CFreal m_Cp;

  /// speed of sound squared
  CFreal m_a2;

  /// Mach number
  CFreal m_Mach;

  /// Temperature
  CFreal m_T;

  /// Temporary Storage for evaluation of Alpha
  RealVector m_eval;

  /// Storage for lift coeficient on trs
  CFreal m_lift;

  /// Storage for lift coeficient on trs
  CFreal m_drag;

  /// force coefficient in XX
  CFreal m_xForceCoef;

  /// force coefficient in YY
  CFreal m_yForceCoef;

  /// Incidence of the TRS in degrees
  CFreal m_alphadeg;

  /// Incidence of the TRS in radians
  CFreal m_alpharad;

  /// velocity at infinity
  CFreal m_uInf;

  /// density at infinity
  CFreal m_rhoInf;

  /// pressure at infinity
  CFreal m_pInf;

  ///flag for appending iteration
  bool m_appendIter;

  ///flag for appending time
  bool m_appendTime;

}; /// end of class Euler2DConsComputeAeroFVMCC

//////////////////////////////////////////////////////////////////////////////

    } /// namespace AeroCoef

} /// namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif /// COOLFluiD_Numerics_AeroCoef_Euler2DConsComputeAeroFVMCC_hh
