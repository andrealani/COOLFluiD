#ifndef COOLFluiD_Physics_Heat_Heat2DDiffusivePrim_hh
#define COOLFluiD_Physics_Heat_Heat2DDiffusivePrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "Heat2DDiffusiveVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Heat {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a Heat physical model 2D for primitive
   * variables
   *
   * @author Andrea Lani
   */
class Heat2DDiffusivePrim : public Heat2DDiffusiveVarSet {
public:

  /**
   * Constructor
   */
  Heat2DDiffusivePrim(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model);

  /**
   * Default destructor
   */
  ~Heat2DDiffusivePrim();

}; // end of class Heat2DDiffusivePrim

//////////////////////////////////////////////////////////////////////////////

    } // namespace Heat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_Heat_Heat2DDiffusive_hh
