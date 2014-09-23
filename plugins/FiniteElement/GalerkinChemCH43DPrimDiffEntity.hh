#ifndef COOLFluiD_Numerics_FiniteElement_GalerkinChemCH43DPrimDiffEntity_hh
#define COOLFluiD_Numerics_FiniteElement_GalerkinChemCH43DPrimDiffEntity_hh

//////////////////////////////////////////////////////////////////////////////

#include "DiffusiveEntity.hh"
#include "Framework/State.hh"
#include "Framework/GeometricEntity.hh"
#include "MathTools/RealMatrix.hh"
#include "Chemistry/CH4/ChemCH43DPrimDiffusive.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class represents a IntegrableEntity for a Galerkin Diffusive Term
   *
   * @author Tiago Quintino
   */
class GalerkinChemCH43DPrimDiffEntity : public DiffusiveEntity {

public:

  /// Default constructor without arguments
  GalerkinChemCH43DPrimDiffEntity(const std::string& name);

  /// Default destructor
  ~GalerkinChemCH43DPrimDiffEntity();

  /**
   * Setup the strategy
   */
  virtual void setup();

  /**
   * Overloading of operator()
   */
  virtual RealMatrix& operator()();

private:

  /// Variable Set
  Common::SafePtr<COOLFluiD::Physics::Chemistry::CH4::ChemCH43DPrimDiffusive> _chemDiffVarSet;

}; // end of class GalerkinChemCH43DPrimDiffEntity

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement
  } // namespace Numerics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteElement_GalerkinChemCH43DPrimDiffEntity_hh
