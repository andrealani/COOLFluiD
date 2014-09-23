#ifndef COOLFluiD_UFEM_UFEMTerm_hh
#define COOLFluiD_UFEM_UFEMTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ConcreteProvider.hh"

#include "UFEM/UFEMSolverData.hh"
#include "UFEM/AssemblyData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace UFEM {

  class ElemProps;

//////////////////////////////////////////////////////////////////////////////

class UFEM_API UFEMTerm : public UFEMSolverStrategy {

public: // typedefs

  /// number of contructor arguments
  enum { NARGS = 2 };
  /// provider for UFEMTerm's
  typedef Environment::ConcreteProvider < UFEMTerm, UFEMTerm::NARGS > PROVIDER;
  /// first argument is name of object
  typedef const std::string& ARG1;
  /// second argument is name the ElemProps
  typedef Common::SafePtr<ElemProps> ARG2;

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Gets the Class name
  static std::string getClassName() { return "UFEMTerm"; }

  /// Constructor.
  explicit UFEMTerm (const std::string& name, Common::SafePtr<ElemProps> props );

  /// Destructor.
  virtual ~UFEMTerm ();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure (Config::ConfigArgs& args);

  /// Set up private data and data of the aggregated classes
  virtual void setup ();

  /// Deallocate the private data
  virtual void unsetup ();

  /// preprocessing (clearing some counters, preparations, etc)
  virtual void pre () {};

  /// Computes the contribution of this term and adds it to the
  /// T; A; b; dTdu.u; dAdu.u; dbdu
  virtual void compute ( const Framework::GeometricEntity& cell , AssemblyData& adata ) = 0;

  /// postprocessing (printing some features, etc)
  virtual void post () {};

}; // class UFEMTerm

//////////////////////////////////////////////////////////////////////////////

  } // namespace UFEM
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

UFEM_Factory(UFEMTerm) // define the factory instance

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_UFEMTerm_hh

