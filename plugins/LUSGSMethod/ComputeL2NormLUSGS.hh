#ifndef COOLFluiD_LUSGSMethod_ComputeL2NormLUSGS_hh
#define COOLFluiD_LUSGSMethod_ComputeL2NormLUSGS_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/GlobalReduce.hh"

#include "LUSGSMethod/ComputeNormLUSGS.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace LUSGSMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class computes the L2 norm of the residuals (solution updates) when using the LUSGS method.
 *
 * @author Kris Van den Abeele
 * @author Matteo Parsani
 *
 */
class ComputeL2NormLUSGS : public ComputeNormLUSGS {

public:

  /// For global reduce
  typedef CFreal GR_RESULTTYPE;

  /// Default constructor without arguments
  ComputeL2NormLUSGS(const std::string& name);

  /// Default destructor
  virtual ~ComputeL2NormLUSGS();

  /// Calculates the norms
  RealVector compute ();

  /**
   * Adds contribution of the current states set to the residuals.
   */
  void addStatesSetContribution();

  /// Retrieves the value for the global reduce of the result
  CFreal GR_GetLocalValue () const;

  /// Global reduce of the result
  static void GR_Combine (const CFreal & S1, const CFreal & S2, CFreal & S3);

  /// Setup the object
  virtual void setup();

  /**
   * Returns the DataSocket's that this numerical strategy needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

protected: // data

  /// object to perform the global reduce
  Framework::GlobalReduce<ComputeL2NormLUSGS> m_gr;
  
  /// socket for rhs of current set of states
  Framework::DataSocketSink< CFreal > socket_rhsCurrStatesSet;

}; // end of class ComputeL2NormLUSGS

//////////////////////////////////////////////////////////////////////////////

    } // namespace LUSGSMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_LUSGSMethod_ComputeL2NormLUSGS_hh
