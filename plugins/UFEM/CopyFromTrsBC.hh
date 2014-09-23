#ifndef COOLFluiD_UFEM_CopyFromTrsBC_hh
#define COOLFluiD_UFEM_CopyFromTrsBC_hh

//////////////////////////////////////////////////////////////////////////////

#include "UFEM/UFEMSolverData.hh"
#include "UFEM/DirichletBC.hh"

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
 * This class represents a CopyFromTrs Boundary condition command
 *
 * @author Tamas Banyai
 *
 */
class UFEM_API CopyFromTrsBC : public DirichletBC {
public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor
  CopyFromTrsBC(const std::string& name);

  /// Destructor
  ~CopyFromTrsBC() {}

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  virtual void setup();

  /// Unsetup the private data and data of the aggregated classes
  /// in this command after the processing phase
  virtual void unsetup() {}

  /// Configures the command
  virtual void configure ( Config::ConfigArgs& args );

protected:

  /// Execute on the current TRS
  virtual void executeOnTrs();

  /// Calculating the variable values
  void computeStateValuesCopyFromTrsBC(const Framework::State* currState);

private:

  /// If nonzero, it looks up the trs specified and the values taken from that trs are being set.
  /// the Def field at the BC in CFcase file will act as a controller, it should be FI for just taking where FI is the name of the variable.
  /// and for example 1.01*FI takes 1.01 times the var from the source trs.
  std::string m_copyFromTrs;

  /// This flag sets whether the states of the source (either from file or from trs) matches in number and in order with current trs.
  /// @todo: it is hardcoded to true currently, implement interpolator (OldLookupTable and OldLookupTable2D -> lookup table for 1D,2D,3D needed)
  bool m_copyMatches;

}; // end of class CopyFromTrsBC

//////////////////////////////////////////////////////////////////////////////

  } // namespace UFEM

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_CopyFromTrsBC_hh
