#ifndef COOLFluiD_UFEM_ElemProps_hh
#define COOLFluiD_UFEM_ElemProps_hh

//////////////////////////////////////////////////////////////////////////////

#include "Environment/ConcreteProvider.hh"

#include "UFEM/UFEMSolverData.hh"

//////////////////////////////////////////////////////////////////////////////
namespace COOLFluiD {
    namespace UFEM {

//////////////////////////////////////////////////////////////////////////////

class UFEM_API ElemProps :  public UFEMSolverStrategy {

public: // typedefs

  /// number of contructor arguments
  enum { NARGS = 1 };
  /// provider for UFEMTerm's
  typedef Environment::ConcreteProvider < ElemProps, ElemProps::NARGS > PROVIDER;
  /// first argument is name of object
  typedef const std::string& ARG1;

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Gets the Class name
  static std::string getClassName() { return "ElemProps"; }

  /// Constructor.
  explicit ElemProps (const std::string& name);

  /// Destructor.
  virtual ~ElemProps ();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure (Config::ConfigArgs& args);

  /// Set up private data and data of the aggregated classes
  virtual void setup ();

  /// Deallocate the private data
  virtual void unsetup ();

  /// Prepare the computation by calculating a priori some
  /// element properties like volume, normals, etc...
  virtual void prepare ( const Framework::GeometricEntity& geoent ) = 0;

}; // class ElemProps

//////////////////////////////////////////////////////////////////////////////

    } // namespace UFEM
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

UFEM_Factory(ElemProps) // define the factory

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_ElemProps_hh

