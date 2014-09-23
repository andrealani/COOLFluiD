#ifndef COOLFluiD_UFEM_QuadP1P1Cell_HeatConduction_hh
#define COOLFluiD_UFEM_QuadP1P1Cell_HeatConduction_hh

//////////////////////////////////////////////////////////////////////////////

#include "UFEM/UFEMTerm.hh"
#include "UFEM/AssemblyData.hh"

#include "UFEM/QuadP1P1Cell/UFEMQuadP1P1Cell.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace UFEM {
    namespace QuadP1P1Cell {

//////////////////////////////////////////////////////////////////////////////

class UFEMQuadP1P1Cell_API HeatConduction : public UFEMTerm {

public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Gets the Class name
  static std::string getClassName() { return "HeatConduction"; }

  /// Constructor.
  explicit HeatConduction ( const std::string& name, Common::SafePtr<ElemProps> props );

  /// Destructor.
  virtual ~HeatConduction ();

  /// Returns the DataSocket's that this command needs as sinks
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();

  /// Returns the DataSocket's that this command provides as sources
  /// @return a vector of SafePtr with the DataSockets
  virtual std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure (Config::ConfigArgs& args);

  /// Set up private data and data of the aggregated classes
  virtual void setup ();

  /// Deallocate the private data
  virtual void unsetup ();

  /// Computes the contribution of this term and adds it to the
  /// T; A; b; dTdu.u; dAdu.u; dbdu
  virtual void compute ( const Framework::GeometricEntity& cell , AssemblyData& adata );

private: // data

  /// Heat Conductivity
  CFreal m_a;

  /// Use of the framework integration
  bool m_cf_integrator;

protected:

  /// socket for states
  Framework::DataSocketSink< Framework::State* > socket_interStates;


}; // class HeatConduction

//////////////////////////////////////////////////////////////////////////////

    }  // end namespace QuadP1P1Cell
  }  // end namespace UFEM
}  // end namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_QuadP1P1Cell_HeatConduction_hh

