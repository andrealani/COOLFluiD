#ifndef COOLFluiD_Numerics_FiniteElement_CoupledNeumannEntity_hh
#define COOLFluiD_Numerics_FiniteElement_CoupledNeumannEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "NeumannEntity.hh"
#include "Framework/Storage.hh"
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
class CoupledNeumannEntity : public NeumannEntity {
public:

  /// Type for the provider of this abstract class
  typedef Framework::BaseMethodStrategyProvider<FiniteElementMethodData,CoupledNeumannEntity> PROVIDER;


 /**
  * Default constructor without arguments
  * @see COOLFluiD::Framework::IntegrableEntity()
  */
  CoupledNeumannEntity(const std::string& name);

 /**
  * Default destructor
  */
  ~CoupledNeumannEntity();

 /**
  * Set up the member data
  */
  void setup();

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
  Common::SafePtr<Framework::VectorialFunction> getRobinVectorialFunction() const
  {
    cf_assert(_vfRobin.isNotNull());
    return _vfRobin;
  }

  /**
   * Overloading of operator()
   * @pre must have first setVarSet()
   */
  virtual RealVector& operator()() = 0;


protected: // data

  /// DataHandle index
  CFuint _index;

  /// IsAccepted DataHandle index
  CFuint _isAcceptedIdx;

  /// DataHandle containing the values at the gauss points
  Framework::DataHandle< RealVector> _interfaceData;

  /// DataHandle containing the flag of acceptance of the gauss points
  Framework::DataHandle< CFreal> _isAccepted;

  /// vector storing the variables to pass to the vectorial function
  RealVector _varsRobin;

  /// acquaintance of the VectorialFunction for the Neumann entity
  Common::SafePtr<Framework::VectorialFunction> _vfRobin;

  /// vector storing the temporary face normal
  bool _isRobin;

}; // end of class CoupledNeumannEntity

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteElement

} // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_CoupledNeumannEntity_hh
