#ifndef COOLFluiD_UFEM_TriagP1P1Cell_CellProps_hh
#define COOLFluiD_UFEM_TriagP1P1Cell_CellProps_hh

//////////////////////////////////////////////////////////////////////////////

#include "UFEM/ElemProps.hh"
#include "UFEM/TriagP1P1Cell/UFEMTriagP1P1Cell.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace UFEM {
    namespace TriagP1P1Cell {

//////////////////////////////////////////////////////////////////////////////

class UFEMTriagP1P1Cell_API CellProps : public ElemProps {

public: // functions

  struct UFEMTriagP1P1Cell_API CellData
  {
    /// Volume of the cell
    CFreal vol;
    /// X components of the inward normal vectors of the faces
    /// Note: the magnitude of the normal vectors are equivalent with the size of that face
    CFreal nx[3];
    /// Y components of the inward normal vectors of the faces
    /// Note: the magnitude of the normal vectors are equivalent with the size of that face
    CFreal ny[3];
  };

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

  /// Returns constant reference to  TriagP1P1Cell cell data
  const CellData& getCellData() const { return m_celldata; }

private: // data

  /// TriagP1P1Cell cell data
  CellData m_celldata;

}; // class CellProps

//////////////////////////////////////////////////////////////////////////////

    }  // end namespace TriagP1P1Cell
  }  // end namespace UFEM
}  // end namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_UFEM_TriagP1P1Cell_CellProps_hh

