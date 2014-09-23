#include "MHD/MHD.hh"
#include "MHD2DProjectionLinearPrim.hh"
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

Environment::ObjectProvider<MHD2DProjectionLinearPrim, JacobianLinearizer, MHDModule, 1> mhd2DProjectionLinearPrimProvider("MHD2DProjectionLinearPrim");

//////////////////////////////////////////////////////////////////////////////

MHD2DProjectionLinearPrim::MHD2DProjectionLinearPrim(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
  _model(model->getImplementor()->getConvectiveTerm().d_castTo<MHDProjectionTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

MHD2DProjectionLinearPrim::~MHD2DProjectionLinearPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DProjectionLinearPrim::linearize(const std::vector<State*>&
					  statesInCell)
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

  linearData[MHDProjectionTerm::RHO] = _avZ[0];
  linearData[MHDProjectionTerm::VX] = _avZ[1];
  linearData[MHDProjectionTerm::VY] = _avZ[2];
  linearData[MHDProjectionTerm::VZ] = _avZ[3];
  linearData[MHDProjectionTerm::V] = sqrt(linearData[MHDProjectionTerm::VX]*
					  linearData[MHDProjectionTerm::VX] +
					  linearData[MHDProjectionTerm::VY]*
					  linearData[MHDProjectionTerm::VY] +
					  linearData[MHDProjectionTerm::VZ]*
					  linearData[MHDProjectionTerm::VZ]);

  linearData[MHDProjectionTerm::P] = _avZ[7];
  linearData[MHDProjectionTerm::A] = sqrt(_model->getGamma()*
					  linearData[MHDProjectionTerm::P]/
					  linearData[MHDProjectionTerm::RHO]);
  linearData[MHDProjectionTerm::BX] = _avZ[4];
  linearData[MHDProjectionTerm::BY] = _avZ[5];
  linearData[MHDProjectionTerm::BZ] = _avZ[6];
  linearData[MHDProjectionTerm::B] = sqrt(linearData[MHDProjectionTerm::BX]*
					  linearData[MHDProjectionTerm::BX] +
					  linearData[MHDProjectionTerm::BY]*
					  linearData[MHDProjectionTerm::BY] +
					  linearData[MHDProjectionTerm::BZ]*
					  linearData[MHDProjectionTerm::BZ]);
  linearData[MHDProjectionTerm::PHI] = _avZ[8];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
