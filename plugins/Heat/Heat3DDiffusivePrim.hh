#ifndef COOLFluiD_Physics_Heat_Heat3DDiffusivePrim_hh
#define COOLFluiD_Physics_Heat_Heat3DDiffusivePrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "Heat3DDiffusiveVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Heat {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Heat physical model 3D for primitive
   * variables
   *
   * @author Andrea Lani
   */
class Heat3DDiffusivePrim : public Heat3DDiffusiveVarSet {
public:

  /**
   * Constructor
   */
  Heat3DDiffusivePrim(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Heat3DDiffusivePrim();

}; // end of class Heat3DDiffusivePrim

//////////////////////////////////////////////////////////////////////////////

    } // namespace Heat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Heat_Heat3DDiffusive_hh
