#ifndef COOLFluiD_Numerics_FluctSplit_StrongSymmetryPlaneImpl_hh
#define COOLFluiD_Numerics_FluctSplit_StrongSymmetryPlaneImpl_hh

//////////////////////////////////////////////////////////////////////////////

#include "StrongImplBC.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/// This class represents a strong slip wall bc for Euler2D
/// @author Andrea Lani
class FluctSplit_API StrongSymmetryPlaneImpl : public StrongImplBC {

public:

  /// Defines the Config Option's of this class
  /// @param options a OptionList where to add the Option's
  static void defineConfigOptions(Config::OptionList& options);

  /// Constructor.
  StrongSymmetryPlaneImpl(const std::string& name);

  /// Default destructor
  ~StrongSymmetryPlaneImpl();

  /// Set up private data and data of the aggregated classes
  /// in this command before processing phase
  void setup();

protected:

  /// Execute on a set of dofs
  void executeOnTrs();

private:

  /// array specifying the components to annull
  std::vector<CFuint> _varIDToAnnull;

}; // end of class StrongSymmetryPlaneImpl

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_StrongSymmetryPlaneImpl_hh
