#include "GalerkinChemCH43DPrimDiffEntity.hh"
#include "DiffusiveEntity.hh"
#include "Environment/ObjectProvider.hh"
#include "MathTools/RealMatrix.hh"
#include "FiniteElement/FiniteElement.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::Chemistry::CH4;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace Numerics {
    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<GalerkinChemCH43DPrimDiffEntity, FiniteElementMethodData,DiffusiveEntity, FiniteElementModule> GalerkinChemCH43DPrimDiffEntityProvider("GalerkinChemCH43DDiffusivePrim");

//////////////////////////////////////////////////////////////////////////////

GalerkinChemCH43DPrimDiffEntity::GalerkinChemCH43DPrimDiffEntity(const std::string& name) :
  DiffusiveEntity(name),
  _chemDiffVarSet(CFNULL)
{
}

//////////////////////////////////////////////////////////////////////////////

GalerkinChemCH43DPrimDiffEntity::~GalerkinChemCH43DPrimDiffEntity()
{
}

//////////////////////////////////////////////////////////////////////////////

void GalerkinChemCH43DPrimDiffEntity::setup()
{
  DiffusiveEntity::setup();

  _chemDiffVarSet = _diffVarSet.d_castTo<ChemCH43DPrimDiffusive>();

}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& GalerkinChemCH43DPrimDiffEntity::operator() ()
{

  const CFuint iState = _localElemData->iState;
  const CFuint jState = _localElemData->jState;

  const CFuint iQuadPoint =  _localElemData->quadPointID;
  const RealMatrix& grad = (*_localElemData->gradValues)[iQuadPoint];

  std::vector<CFreal>& vecD = _chemDiffVarSet->getDiffusiveCoefs();

  for(CFuint i = 0; i < vecD.size(); ++i) {
    _result(i,i) = vecD[i] * (grad(iState,XX)*grad(jState,XX) +
                              grad(iState,YY)*grad(jState,YY) +
                              grad(iState,ZZ)*grad(jState,ZZ));
  }

  return _result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace COOLFluiD
  } // namespace Numerics
} // namespace FiniteElement

