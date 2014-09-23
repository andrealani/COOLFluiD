#ifndef COOLFluiD_Framework_AdaptSubSystem_hh
#define COOLFluiD_Framework_AdaptSubSystem_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/StandardSubSystem.hh"
#include "Framework/MaxTimeCondition.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

    class Method;
    class ConvergenceMethod;
    class SpaceMethod;
    class ErrorEstimatorMethod;
    class CouplerMethod;
    class MeshAdapterMethod;
    class OutputFormatter;
    class DataProcessing;
    class MeshCreator;
    class PhysicalModelImpl;
    class NumericalCommand;
    class NumericalStrategy;
    class ComputeNorm;
    class ComputeCFL;
    class LinearSystemSolver;
    class MeshDataAdapter;

}

  namespace Numerics {

    namespace RemeshMeandros {

//////////////////////////////////////////////////////////////////////////////

/// A AdaptSubSystem is a concrete implementation of the SubSystem.
/// @author Jurek Majewski
class AdaptSubSystem : public Framework::StandardSubSystem
{
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// The default constructor without arguments.
  /// @see ~AdaptSubSystem()
  AdaptSubSystem(const std::string& name);

  /// Destructor
  /// @see AdaptSubSystem()
  ~AdaptSubSystem();

  /// Configures this Simualtion.
  /// Sets up the data for this Object.
  virtual void configure ( Config::ConfigArgs& args );

  /// Run (Process) Phase. All the big number crunching work goes here.
  void run();

  /// Adds the ActionListener's of this EventListener to the EventHandler
  void registActionListeners();

private:

  /// Number of error estimations during solution
  CFuint  m_estimateNum;

  Common::SafePtr<Framework::MaxTimeCondition> m_pMaxTimeStopCond;

}; // class AdaptSubSystem

//////////////////////////////////////////////////////////////////////////////

    } // namespace RemeshMeandros

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_AdaptSubSystem_hh
