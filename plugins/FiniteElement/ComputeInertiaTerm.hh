#ifndef COOLFluiD_Numerics_FiniteElement_ComputeInertiaTerm_hh
#define COOLFluiD_Numerics_FiniteElement_ComputeInertiaTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NullableObject.hh"
#include "Framework/ComputeTerm.hh"
#include "Framework/GeometricEntity.hh"
#include "InertiaEntity.hh"
#include "Common/SafePtr.hh"
#include "FiniteElementMethodData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    class Element;
  }

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents an object computing a Inertia term
 * with Finite Element method
 *
 * @author Thomas Wuilbaut
 * @author Tiago Quintino
 */
class ComputeInertiaTerm : public Framework::ComputeTerm<FiniteElementMethodData>{

public://types

  typedef Framework::BaseMethodStrategyProvider<FiniteElementMethodData,ComputeInertiaTerm> PROVIDER;

public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Constructor
   */
  ComputeInertiaTerm(const std::string& name);

  /**
   * Default destructor
   */
  ~ComputeInertiaTerm();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * This ComputeTerm is null if the InertiaVarSet is null
   */
  bool isNull() const
  {
    return false;
  }

  /**
   * Set up private data to prepare the simulation
   */
  void setup();

  /**
   * Computes the Inertia Term
   */
  void computeTerm(Framework::GeometricEntity* const cell, RealMatrix& result);

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "ComputeInertiaTerm";
  }


protected:

  /// Inertia Integrable Entity
  Common::SafePtr<InertiaEntity> _inertiaEntity;

}; // end of class ComputeInertiaTerm

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_ComputeInertiaTerm_hh
