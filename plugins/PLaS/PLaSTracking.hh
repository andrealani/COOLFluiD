#ifndef COOLFluiD_PLaS_PLaSTracking_hh
#define COOLFluiD_PLaS_PLaSTracking_hh

#include "Framework/DataProcessingMethod.hh"
#include "Framework/SpaceMethod.hh"
#include "PLaS/PLaSTrackingData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace PLaS {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a DataProcessingMethod interfacing the PLaS library
class PLaSTracking : public Framework::DataProcessingMethod {

 public:  // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Default constructor
  PLaSTracking(const std::string& name);

  /// Default destructor
  ~PLaSTracking();

  /// Configure
  void configure(Config::ConfigArgs& args);

  /// Gets the Class name
  static std::string getClassName() { return "PLaSTracking"; }


 protected:  // abstract interface implementations

  /// Execute the data processing
  void processDataImpl() {
    m_process->execute();
  }

  /// UnSets the data of the method
  /// @see Method::setMethod()
  virtual void setMethodImpl() {
    DataProcessingMethod::setMethodImpl();
    setupCommandsAndStrategies();
    m_setup->execute();
  }

  /// Sets up the data for the method commands to be applied
  /// @see Method::unsetMethod()
  virtual void unsetMethodImpl() {
    unsetupCommandsAndStrategies();
    m_unsetup->execute();
    DataProcessingMethod::unsetMethodImpl();
  }

  /// Gets the Data aggregator of this method
  /// @return SafePtr to the MethodData
  Common::SafePtr< Framework::MethodData > getMethodData() const {
    return m_data.getPtr();
  }


 public:  // helper functions

  /// Get PLaS data of the dispersed phase (per node)
  const PLAS_PHASE_DATA* getPhaseData() const {
    return m_data->getPhaseData();
  }

  /// Get PLaS input parameters
  const PLAS_INPUT_PARAM& getInputParameters() const {
    return m_data->getInputParameters();
  }


 public:  // non-standard helper functions

  /// Execute the data processing (public version)
  void process() {
    processDataImpl();
  }

  /// Block further processing
  void block() {
    m_data->m_block = true;
  }

  /// Unblock further processing
  void unblock() {
    m_data->m_block = false;
  }

  /// Set primary phase properties
  /// @param _rho Primary (continuous) phase density
  /// @param _nu Primary (continuous) phase kinematic viscosity
  /// @param _dt Eulerian time scale
  /// @param _cpCont Specific heat coefficient of the flow medium
  /// @param _kCont Thermal conductivity of the flow medium
  /// @param _numunk Number of unknown variables for the primary phase flow
  void setFlowProperties(double _rho, double _nu, double _dt, double _cpCont, double _kCont, int _numunk) {
    m_data->setFlowProperties(_rho,_nu,_dt,_cpCont,_kCont,_numunk);
  }

  /// Set iteration properties
  /// @param _iter Current iteration
  /// @param _time Current time
  /// @param _dt Iteration time-step
  /// @param _output Flag for writing output
  /// @param _numExtEnt Number of bubbles coming from external code
  /// @param _extEntPos Positions of bubbles coming from external code
  /// @param _extEntVel Velocities of bubbles coming from external code
  /// @param _extEntTemp Temperature of bubbles coming from external code
  /// @param _extEntDiam Diameters of bubbles coming from external code
  void setIterationProperties(int _iter, double _time, double _dt, int _output, int _numExtEnt=0, double *_extEntPos=CFNULL, double *_extEntVel=CFNULL, double *_extEntTemp=CFNULL, double *_extEntDiam=CFNULL) {
    m_data->setIterationProperties(_iter, _time, _dt, _output, _numExtEnt, _extEntPos, _extEntVel, _extEntTemp, _extEntDiam);
  }


 private:  // data

  /// String for the setup command
  std::string m_setupStr;

  /// String for the unsetup command
  std::string m_unsetupStr;

  /// Strings for the process command
  std::string m_processStr;

  /// Pointer to the setup command
  Common::SelfRegistPtr< PLaSTrackingCom > m_setup;

  /// Pointer to the unsetup command
  Common::SelfRegistPtr< PLaSTrackingCom > m_unsetup;

  /// Pointer to the process command
  Common::SelfRegistPtr< PLaSTrackingCom > m_process;

  /// The data to share between PLaSTracking commands
  Common::SharedPtr< PLaSTrackingData > m_data;

};

//////////////////////////////////////////////////////////////////////////////

  }  // namespace PLaS
}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_PLaS_PLaSTracking_hh

