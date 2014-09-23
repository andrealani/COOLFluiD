#ifndef COOLFluiD_Numerics_FiniteElement_CoupledNeumannBCEntity_hh
#define COOLFluiD_Numerics_FiniteElement_CoupledNeumannBCEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/IntegrableEntity.hh"
#include "Framework/State.hh"
#include "MathTools/RealMatrix.hh"
#include "Common/SafePtr.hh"
#include "Framework/Storage.hh"
#include "Framework/VectorialFunction.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a IntegrableEntity for a CoupledNeumannBC
   *
   * @author Thomas Wuilbaut
   *
   * @see IntegrableEntity
   *
   */
class CoupledNeumannBCEntity : public COOLFluiD::Framework::IntegrableEntity {
public:

  typedef Environment::ConcreteProvider<CoupledNeumannBCEntity> PROVIDER;

 /**
  * Default constructor without arguments
  * @see COOLFluiD::Framework::IntegrableEntity()
  */
  CoupledNeumannBCEntity();

 /**
  * Default destructor
  */
  ~CoupledNeumannBCEntity();

 /**
  * Set up the member data
  */
  void setup();

 /**
  * Set the state index
  */
  void setStateIndex(const CFuint i)
  {
    _iState = i;
  }

 /**
  * Set the dataHandle index
  */
  void setDataHandleIndex(const CFuint idx)
  {
    _index = idx;
  }

 /**
  * Set the dataHandle index
  */
  void setIsAcceptedDataHandleIndex(const CFuint isAcceptedIdx)
  {
    _isAcceptedIdx = isAcceptedIdx;
  }

  /**
  * Set the dataHandle containing the values at the gauss points
  */
  void setDataHandle(Framework::DataHandle< RealVector> interfaceData)
  {
    _interfaceData = interfaceData;
  }


 /**
  * Set the dataHandle containing the flag for accepted values
  */
  void setIsAcceptedDataHandle(Framework::DataHandle< CFreal> isAccepted)
  {
    _isAccepted = isAccepted;
  }

 /**
  * Set the vectorial function
  */
  void setVectorialFunction(Common::SafePtr<Framework::VectorialFunction> vf)
  {
    _vf = vf;
    cf_assert(_vf.isNotNull());
  }

 /**
  * Set the vectorial function for RobinBC
  */
  void setVectorialFunctionRobin(Common::SafePtr<Framework::VectorialFunction> vf)
  {
    _vfRobin = vf;
    cf_assert(_vfRobin.isNotNull());
  }

 /**
  * Set the vectorial function for RobinBC
  */
  void setIsRobinBC(bool isRobin, const CFuint size)
  {
    _isRobin = isRobin;
    _varsRobin.resize(size);
  }


 /**
  * Get the vectorial function
  */
  Common::SafePtr<Framework::VectorialFunction> getVectorialFunction() const
  {
    cf_assert(_vf.isNotNull());
    return _vf;
  }

 /**
  * Get the vectorial function
  */
  Common::SafePtr<Framework::VectorialFunction> getRobinVectorialFunction() const
  {
    cf_assert(_vfRobin.isNotNull());
    return _vfRobin;
  }

  /**
   * Overloading of operator()
   */
  RealVector& operator()(const Framework::State& state, const RealVector& shapeFW, Framework::GeometricEntity* geo);

  /**
   * Get the size of the result
   */
  CFuint size() { return _result.size(); }

protected: // data

  /// vector storing the temporary computation
  RealVector _result;

  /// DataHandle index
  CFuint _index;

  /// IsAccepted DataHandle index
  CFuint _isAcceptedIdx;

  /// State index
  CFuint _iState;

  /// DataHandle containing the values at the gauss points
  Framework::DataHandle< RealVector> _interfaceData;

  /// DataHandle containing the flag of acceptance of the gauss points
  Framework::DataHandle< CFreal> _isAccepted;

  /// dimension
  CFuint _dim;

  /// nb equations
  CFuint _nbEq;

  /// vector storing the variables to pass to the vectorial function
  RealVector _vars;

  /// vector storing the variables to pass to the vectorial function
  RealVector _varsRobin;

  /// acquaintance of the VectorialFunction for the Neumann entity
  Common::SafePtr<Framework::VectorialFunction> _vf;

  /// acquaintance of the VectorialFunction for the Neumann entity
  Common::SafePtr<Framework::VectorialFunction> _vfRobin;

  /// vector storing the temporary face normal
  RealVector _normal;

  /// vector storing the temporary face normal
  bool _isRobin;

}; // end of class CoupledNeumannBCEntity

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteElement

} // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_CoupledNeumannBCEntity_hh
