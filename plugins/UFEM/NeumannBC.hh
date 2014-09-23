#ifndef COOLFluiD_UFEM_NeumannBC_hh
#define COOLFluiD_UFEM_NeumannBC_hh

//////////////////////////////////////////////////////////////////////////////

#include "UFEM/UFEMSolverData.hh"
#include "UFEM/BaseBC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework
  {
    template< typename, typename > class DataSocketSink;
    class VectorialFunction;
  }

  namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a ImplicitNeumann Boundary condition command
 *
 * @author Tiago Quintino
 * @author reworked for UFEM by Tamas Banyai
 *
 */
class UFEM_API NeumannBC : public BaseBC {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  NeumannBC(const std::string& name);

  /// Destructor
  ~NeumannBC() {}

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Unsetup the private data and data of the aggregated classes
  /// in this command after the processing phase
  virtual void unsetup() {}

  /// Configures the command
  virtual void configure ( Config::ConfigArgs& args );

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector< Common::SafePtr< Framework::BaseDataSocketSink > > needsSockets();

protected:

  /// Execute on the current TRS
  virtual void executeOnTrs();

  /// Calculating the variable values
  void computeStateValuesNeumannBC(const Framework::State* currState);

private:

  /// Method is Implicit or Explicit (adding to matrix or rhs).
  /// If implicit is false: dU/dn*k_U=Vars_U where Vars_U is the value specified in the BC's Vars field and k=diffusion coefficient of the equation (remember the weak formulation).
  /// If implicit is true: dU/dn*k_U=U*Vars_U (so Vars_U multiplied by U is the equivalent with implicit=false).
  bool m_isImplicit;

  /// Socket for RHS
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// Socket for isUpdated
  Framework::DataSocketSink<bool> socket_isUpdated;

  /// Socket for neighbor states
  Framework::DataSocketSink< std::valarray< Framework::State* > > socket_bStatesNeighbors;

}; // end of class NeumannBC

//////////////////////////////////////////////////////////////////////////////

  } // namespace UFEM

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_NeumannBC_hh
