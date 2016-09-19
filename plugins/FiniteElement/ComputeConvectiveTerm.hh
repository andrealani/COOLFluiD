#ifndef COOLFluiD_Numerics_FiniteElement_ComputeConvectiveTerm_hh
#define COOLFluiD_Numerics_FiniteElement_ComputeConvectiveTerm_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NullableObject.hh"
#include "Framework/ComputeTerm.hh"
#include "Framework/GeometricEntity.hh"
#include "ConvectiveEntity.hh"
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
 * This class represents an object computing a convective term
 * with Finite Element method
 */
class ComputeConvectiveTerm :
  public Framework::ComputeTerm<FiniteElementMethodData>{

public://types

  typedef Framework::BaseMethodStrategyProvider<FiniteElementMethodData,ComputeConvectiveTerm> PROVIDER;

public:

  /// Constructor
  ComputeConvectiveTerm(const std::string& name);

  /// Default destructor
  ~ComputeConvectiveTerm();

  /// Configure the object
  virtual void configure ( Config::ConfigArgs& args );

  /// This ComputeTerm is not null
  bool isNull() const
  {
    return false;
  }

  /// Set up private data to prepare the simulation
  void setup(){

    _convectiveEntity = getMethodData().getConvectiveEntity();
    cf_assert(_convectiveEntity.isNotNull());
  }

  /// Computes the Convective Term
  void computeTerm(Framework::GeometricEntity* const cell, RealMatrix& result);

  /// Gets the Class name
  static std::string getClassName()
  {
    return "ComputeConvectiveTerm";
  }
  
  /// Gets the polymorphic type name
  virtual std::string getPolymorphicTypeName() {return getClassName();}
  
protected:

  /// Convective Integrable Entity
  Common::SafePtr<ConvectiveEntity> _convectiveEntity;

}; // end of class ComputeConvectiveTerm


//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_ComputeConvectiveTerm_hh

