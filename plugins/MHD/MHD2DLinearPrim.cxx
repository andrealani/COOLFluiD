#include "MHD/MHD.hh"
#include "MHD2DLinearPrim.hh"
#include "MHDPhysicalModel.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MHD2DLinearPrim, JacobianLinearizer, MHDModule, 1> mhd2DLinearPrimProvider("MHD2DLinearPrim");

//////////////////////////////////////////////////////////////////////////////

MHD2DLinearPrim::MHD2DLinearPrim(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
  _model(model->getImplementor()->getConvectiveTerm().d_castTo<MHDTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

MHD2DLinearPrim::~MHD2DLinearPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DLinearPrim::linearize(const std::vector<State*>& statesInCell)
{
  RealVector& linearData = _model->getPhysicalData();

  // this comes from outside
  // const vector<State*> *const tStates = _transformer.transform
  // (const_cast<std::vector<State*>*>(&statesInCell));

  const CFuint nbStates = statesInCell.size();
  // linearizer does this (numerics dependent)
  _sumZ = 0.;
  for (CFuint iEq = 0; iEq < PhysicalModelStack::getActive()->getNbEq(); ++iEq) {
    for (CFuint iState = 0; iState < nbStates; ++iState) {
      _sumZ[iEq] += (*statesInCell[iState])[iEq];
    }
  }
  _avZ = _sumZ/static_cast<CFreal>(nbStates);
  //////

  linearData[MHDTerm::RHO] = _avZ[0];
  linearData[MHDTerm::VX] = _avZ[1];
  linearData[MHDTerm::VY] = _avZ[2];
  linearData[MHDTerm::VZ] = _avZ[3];
  linearData[MHDTerm::V] = sqrt(linearData[MHDTerm::VX]*linearData[MHDTerm::VX] +
				linearData[MHDTerm::VY]*linearData[MHDTerm::VY] +
				linearData[MHDTerm::VZ]*linearData[MHDTerm::VZ]);

  linearData[MHDTerm::P] = _avZ[7];
  linearData[MHDTerm::A] = sqrt(_model->getGamma()*linearData[MHDTerm::P]/
				linearData[MHDTerm::RHO]);
  linearData[MHDTerm::BX] = _avZ[4];
  linearData[MHDTerm::BY] = _avZ[5];
  linearData[MHDTerm::BZ] = _avZ[6];
  linearData[MHDTerm::B] = sqrt(linearData[MHDTerm::BX]*linearData[MHDTerm::BX] +
				linearData[MHDTerm::BY]*linearData[MHDTerm::BY] +
				linearData[MHDTerm::BZ]*linearData[MHDTerm::BZ]);
}

//////////////////////////////////////////////////////////////////////////////

} // namespace MHD
} // namespace Physics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
