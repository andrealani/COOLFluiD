#ifndef COOLFluiD_Physics_LinEuler_CreateDampingCoefficient_hh
#define COOLFluiD_Physics_LinEuler_CreateDampingCoefficient_hh

//////////////////////////////////////////////////////////////////////////////

//#include "MathTools/FunctionParser.hh"
#include "Framework/DataProcessingData.hh"
//#include "Framework/DataSocketSource.hh"
//#include "Framework/DataSocketSink.hh"
//#include "Framework/VectorialFunction.hh"

//#include "Framework/PhysicalModel.hh"
//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

   namespace LinearizedEuler {

//////////////////////////////////////////////////////////////////////////////

///
/// This class create a socket for the damping coefficient
/// @author Erik Torres


class CreateDampingCoefficient : public Framework::DataProcessingCom {
public: // functions

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor.
   */
  CreateDampingCoefficient(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~CreateDampingCoefficient();

  /** Do I need this?
   * Configure the command
   */
  //virtual void configure (const Config::ConfigArgs& args);

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
  virtual void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command
   */
  virtual void unsetup();

protected: // functions

  /**
   * Execute on a set of dofs
   */
  void executeOnTrs();

private: // data

  /// the socket to the data handle of the state's
  Framework::DataSocketSink < Framework::State* , Framework::GLOBAL > socket_states;

  /// the socket stores the data of the damping coefficient
  Framework::DataSocketSource<CFreal> socket_dampingCoeff;

  /// Maximum damping coefficient (at boundary)
  CFreal m_nuMax;

  /// Radius of influence of the damping zone
  CFreal m_r0;

  /// Exponent for the damping function
  CFreal m_beta;

  ///Name of the boundary associated to the damping zone
  std::vector<std::string> m_BoundaryTRS;

}; // end of class CreateDampingCoefficient

//////////////////////////////////////////////////////////////////////////////

    } // namespace LinEuler

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LinEuler_CreateDampingCoefficient_hh
