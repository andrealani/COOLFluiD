#ifndef COOLFluiD_Numerics_FiniteElement_NeumannEntity_hh
#define COOLFluiD_Numerics_FiniteElement_NeumannEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "FEMIntegrableEntity.hh"
#include "Framework/State.hh"
#include "MathTools/RealMatrix.hh"
#include "Common/SafePtr.hh"
#include "Framework/VectorialFunction.hh"
#include "LocalElementData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a IntegrableEntity for a diffusive term
   *
   * @author Thomas Wuilbaut
   * @author Tiago Quintino
   *
   * @see IntegrableEntity
   *
   */
class NeumannEntity : public FiniteElement::FEMIntegrableEntity {
public:

  /// Type for the provider of this abstract class
  typedef Framework::BaseMethodStrategyProvider<FiniteElementMethodData,NeumannEntity> PROVIDER;

 /**
  * Default constructor without arguments
  * @see COOLFluiD::Framework::IntegrableEntity()
  */
  NeumannEntity(const std::string& name);

 /**
  * Default destructor
  */
  ~NeumannEntity();

 /**
  * Sets the vectorial function
  */
  void setVectorialFunction(Common::SafePtr<Framework::VectorialFunction> vf)
  {
    _vf = vf;
  }

  /**
   * Overloading of operator()
   * @pre must have first setVarSet()
   */
  virtual RealVector& operator()() = 0;

  /// Setup the strategy
  virtual void setup()
  {
    FEMIntegrableEntity::setup();

    _diffVarSet = getMethodData().getDiffusiveVar();
    _nbEqs = Framework::PhysicalModelStack::getActive()->getNbEq();
    _dim = Framework::PhysicalModelStack::getActive()->getDim();
    _result.resize(_nbEqs);
    _normal.resize(_dim);
    _vars.resize(_dim+_dim+_nbEqs+1);
  }

  /**
   * Get the size of the result
   */
  CFuint size() { return _result.size(); }

protected:

  /// Diffusive Variable Set
  Common::SafePtr<Framework::DiffusiveVarSet> _diffVarSet;

  /// matrix storing the temporary computation
  RealVector _result;

  ///  number of equations
  CFuint _nbEqs;

  ///  number of dimension
  CFuint _dim;

  /// acquaintance of the VectorialFunction for the Neumann entity
  Common::SafePtr<Framework::VectorialFunction> _vf;

  /// Vector storing the coordinates, the state values and the time
  RealVector _vars;

  /// Vector storing the normal to the face
  RealVector _normal;

}; // end of class NeumannEntity

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteElement

} // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_NeumannEntity_hh
