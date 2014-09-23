#include "MHD/MHD.hh"
#include "MHD3DProjectionLinearCons.hh"
#include "MHDProjection.hh"
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

Environment::ObjectProvider<MHD3DProjectionLinearCons, JacobianLinearizer, MHDModule, 1> mhd3DProjectionLinearConsProvider("MHD3DProjectionLinearCons");

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionLinearCons::MHD3DProjectionLinearCons(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
  _model(model->getImplementor()->getConvectiveTerm().d_castTo<MHDProjectionTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionLinearCons::~MHD3DProjectionLinearCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionLinearCons::linearize(const std::vector<State*>& statesInCell)
{
  RealVector& linearData = _model->getPhysicalData();

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

  linearData[MHDProjectionTerm::RHO] = _avZ[0];
  linearData[MHDProjectionTerm::VX] = _avZ[1]/_avZ[0];
  linearData[MHDProjectionTerm::VY] = _avZ[2]/_avZ[0];
  linearData[MHDProjectionTerm::VZ] = _avZ[3]/_avZ[0];
  linearData[MHDProjectionTerm::V] = sqrt(linearData[MHDProjectionTerm::VX]*linearData[MHDProjectionTerm::VX] +
                                linearData[MHDProjectionTerm::VY]*linearData[MHDProjectionTerm::VY] +
                                linearData[MHDProjectionTerm::VZ]*linearData[MHDProjectionTerm::VZ]);
  linearData[MHDProjectionTerm::BX] = _avZ[4];
  linearData[MHDProjectionTerm::BY] = _avZ[5];
  linearData[MHDProjectionTerm::BZ] = _avZ[6];
  linearData[MHDProjectionTerm::B] = sqrt(linearData[MHDProjectionTerm::BX]*linearData[MHDProjectionTerm::BX] +
                                linearData[MHDProjectionTerm::BY]*linearData[MHDProjectionTerm::BY] +
                                linearData[MHDProjectionTerm::BZ]*linearData[MHDProjectionTerm::BZ]);

  linearData[MHDProjectionTerm::P] = (_model->getGamma() - 1.0)*
    (_avZ[7] - 0.5*(linearData[MHDProjectionTerm::RHO]*linearData[MHDProjectionTerm::V]*
                    linearData[MHDProjectionTerm::V] + linearData[MHDProjectionTerm::B]*
                    linearData[MHDProjectionTerm::B]));
  linearData[MHDProjectionTerm::A] = sqrt(_model->getGamma()*linearData[MHDProjectionTerm::P]/
                                linearData[MHDProjectionTerm::RHO]);
  linearData[MHDProjectionTerm::PHI] = _avZ[8];
}

//////////////////////////////////////////////////////////////////////////////

} // namespace MHD
} // namespace Physics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
