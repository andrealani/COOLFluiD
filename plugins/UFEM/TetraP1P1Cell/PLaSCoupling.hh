#ifndef COOLFluiD_UFEM_TetraP1P1Cell_PLaSCoupling_hh
#define COOLFluiD_UFEM_TetraP1P1Cell_PLaSCoupling_hh

//////////////////////////////////////////////////////////////////////////////

#include "UFEM/UFEMTerm.hh"
#include "UFEM/AssemblyData.hh"
#include "UFEM/TetraP1P1Cell/CellProps.hh"
#include "PLaS/PLaSTrackingData.hh"


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace UFEM {
    namespace TetraP1P1Cell {

//////////////////////////////////////////////////////////////////////////////

class UFEMTetraP1P1Cell_API PLaSCoupling : public UFEMTerm {

public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Gets the Class name
  static std::string getClassName() { return "PLaSCoupling"; }

  /// Constructor.
  explicit PLaSCoupling ( const std::string& name, Common::SafePtr<ElemProps> props );

  /// Destructor.
  virtual ~PLaSCoupling ();

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
  TetraP1P1Cell::CellProps& m_cellprops;

protected:

	/// pointer to PLaS phase data (pointer to vector per node)
	const PLAS_PHASE_DATA * m_plas_phase_data;

}; // class PLaSCoupling

//////////////////////////////////////////////////////////////////////////////

    }  // end namespace TetraP1P1Cell
  }  // end namespace UFEM
}  // end namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_TetraP1P1Cell_PLaSCoupling_hh

