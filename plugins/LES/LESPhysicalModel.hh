#ifndef COOLFluiD_Physics_LES_LESPhysicalModel_hh
#define COOLFluiD_Physics_LES_LESPhysicalModel_hh

//////////////////////////////////////////////////////////////////////////////

#include "NavierStokes/NavierStokesPhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace LES {

//////////////////////////////////////////////////////////////////////////////

/// This class represents the interface for a LESPhysicalModels
/// @author Kris Van den Abeele
/// @author Ghader Ghorbaniasl
template <int DIM>
 class LESPhysicalModel :
        public Physics::NavierStokes::NavierStokesPhysicalModel<DIM> {
 public:

  /// Constructor without arguments
  LESPhysicalModel(const std::string& name);

  /// Default destructor
  virtual ~LESPhysicalModel();

  /// Configures this object with user parameters
  virtual void configure ( Config::ConfigArgs& args );

  /// Get the name of the Physical Model type
  /// @return name in a std::string
  virtual std::string getTypeName() const
  {
    return std::string("LES" + Common::StringOps::to_str(DIM) + "D");
  }

}; // end of class LESPhysicalModel

//////////////////////////////////////////////////////////////////////////////

    } // namespace LES

} // namespace COOLFluiD

#include "LESPhysicalModel.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_LES_LESPhysicalModel_hh
