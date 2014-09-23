#ifndef COOLFluiD_UFEM_TriagP1P1Cell_MitremBin_hh
#define COOLFluiD_UFEM_TriagP1P1Cell_MitremBin_hh

//////////////////////////////////////////////////////////////////////////////

#include "UFEM/UFEMTerm.hh"
#include "UFEM/AssemblyData.hh"

#include "UFEM/TriagP1P1Cell/CellProps.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace UFEM {
    namespace TriagP1P1Cell {

//////////////////////////////////////////////////////////////////////////////

class UFEMTriagP1P1Cell_API MitremBin : public UFEMTerm {

public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Gets the Class name
  static std::string getClassName() { return "MitremBin"; }

  /// Constructor.
  explicit MitremBin ( const std::string& name, Common::SafePtr<ElemProps> props );

  /// Destructor.
  virtual ~MitremBin ();

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

  /// Cell properties (normals,volume ... etc)
  TriagP1P1Cell::CellProps& m_cellprops;

  /// Charge numbers
	std::vector<CFreal> m_charge;
	
	/// Diffusion constants
  std::vector<CFreal> m_diff;

protected:

  /// socket for states
  Framework::DataSocketSink< Framework::State* > socket_interStates;


}; // class MitremBin

//////////////////////////////////////////////////////////////////////////////

    }  // end namespace TriagP1P1Cell
  }  // end namespace UFEM
}  // end namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_TriagP1P1Cell_MitremBin_hh

