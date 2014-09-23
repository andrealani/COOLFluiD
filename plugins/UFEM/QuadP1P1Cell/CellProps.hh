#ifndef COOLFluiD_UFEM_QuadP1P1Cell_CellProps_hh
#define COOLFluiD_UFEM_QuadP1P1Cell_CellProps_hh

//////////////////////////////////////////////////////////////////////////////

#include "UFEM/ElemProps.hh"
#include "UFEM/QuadP1P1Cell/UFEMQuadP1P1Cell.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace UFEM {
    namespace QuadP1P1Cell {

//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

class UFEMQuadP1P1Cell_API CellProps : public ElemProps {

public: // functions

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Gets the Class name
  static std::string getClassName() { return "CellProps"; }

  /// Constructor.
  explicit CellProps (const std::string& name);

  /// Destructor.
  virtual ~CellProps ();

  /// Configure the data from the supplied arguments.
  /// @param args configuration arguments
  virtual void configure (Config::ConfigArgs& args);

  /// Set up private data and data of the aggregated classes
  virtual void setup ();

  /// Deallocate the private data
  virtual void unsetup ();

  /// Prepare the computation by calculating a priori some
  /// element properties like volume, normals, etc...
  virtual void prepare ( const Framework::GeometricEntity& cell );

}; // class CellProps

//////////////////////////////////////////////////////////////////////////////

    }  // end namespace QuadP1P1Cell
  }  // end namespace UFEM
}  // end namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_QuadP1P1Cell_CellProps_hh

