#ifndef COOLFluiD_UFEM_DirichletBC_hh
#define COOLFluiD_UFEM_DirichletBC_hh

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
 * This class represents a Dirichlet Boundary condition command
 *
 * @author Tiago Quintino
 * @author reworked for UFEM by Tamas Banyai
 *
 */
class UFEM_API DirichletBC : public BaseBC {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  DirichletBC(const std::string& name);

  /// Destructor
  ~DirichletBC() {}

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
  void computeStateValuesDirichletBC(const Framework::State* currState);

  /// Socket for RHS
  Framework::DataSocketSink<CFreal> socket_rhs;

  /// Socket for isUpdated
  Framework::DataSocketSink<bool> socket_isUpdated;

  /// Socket for neighbor states
  Framework::DataSocketSink< std::valarray< Framework::State* > > socket_bStatesNeighbors;

private:

  /// Method is Implicit or Explicit (solving for deltaFI or FI)
  bool m_isImplicit;

  /// Symmetry method to apply
  std::string m_symmetryStr;

  /// ScaleDiagonal symmetry method coefficient
  CFreal m_scale;

}; // end of class DirichletBC

//////////////////////////////////////////////////////////////////////////////

  } // namespace UFEM

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_DirichletBC_hh
