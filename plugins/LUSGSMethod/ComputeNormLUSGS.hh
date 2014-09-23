#ifndef COOLFluiD_LUSGSMethod_ComputeNormLUSGS_hh
#define COOLFluiD_LUSGSMethod_ComputeNormLUSGS_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ComputeNorm.hh"
#include "Framework/DataSocketSink.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a generic interface for computing the norm when using the LUSGS method.
 *
 * @author Kris Van den Abeele
 * @author Matteo Parsani
 *
 */
class ComputeNormLUSGS : public Framework::ComputeNorm {

public:

  typedef Environment::ConcreteProvider<ComputeNormLUSGS,1> PROVIDER;
  typedef const std::string& ARG1;

  /// Default constructor without arguments
  ComputeNormLUSGS(const std::string& name);

  /// Default destructor
  virtual ~ComputeNormLUSGS();

  /// Setup the object
  virtual void setup();

  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /**
   * Sets the residual norms to zero.
   */
  void resetResiduals()
  {
    m_localResiduals = 0.0;
  }

  /**
   * Adds contribution of the current states set to the residuals.
   */
  virtual void addStatesSetContribution() = 0;

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "ComputeNormLUSGS";
  }

protected: // data

  /// socket for current states set index
  Framework::DataSocketSink< CFint > socket_statesSetIdx;

  /// handle to the IDs of the states in each set of states
  Framework::DataSocketSink< std::vector< CFuint > > socket_statesSetStateIDs;

  /// handle to list of booleans telling whether a states set is parallel updatable
  Framework::DataSocketSink< bool > socket_isStatesSetParUpdatable;

  /// local residual norm in this processor
  RealVector m_localResiduals;

  /// number of equations
  CFuint m_nbrEqs;

}; // end of class ComputeNormLUSGS

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_LUSGSMethod_ComputeNormLUSGS_hh
