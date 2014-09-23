#include "MHD/MHD.hh"
#include "MHD2DLinearCons.hh"
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

Environment::ObjectProvider<MHD2DLinearCons, JacobianLinearizer, MHDModule, 1> mhd2DLinearConsProvider("MHD2DLinearCons");

//////////////////////////////////////////////////////////////////////////////

MHD2DLinearCons::MHD2DLinearCons(Common::SafePtr<Framework::PhysicalModel> model) :
  JacobianLinearizer(model),
  _model(model->getImplementor()->
	 getConvectiveTerm().d_castTo<MHDTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

MHD2DLinearCons::~MHD2DLinearCons()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DLinearCons::linearize(const std::vector<State*>& statesInCell)
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

  linearData[MHDTerm::RHO] = _avZ[0];
  linearData[MHDTerm::VX] = _avZ[1]/_avZ[0];
  linearData[MHDTerm::VY] = _avZ[2]/_avZ[0];
  linearData[MHDTerm::VZ] = _avZ[3]/_avZ[0];
  linearData[MHDTerm::V] = sqrt(linearData[MHDTerm::VX]*linearData[MHDTerm::VX] +
				linearData[MHDTerm::VY]*linearData[MHDTerm::VY] +
				linearData[MHDTerm::VZ]*linearData[MHDTerm::VZ]);
  linearData[MHDTerm::BX] = _avZ[4];
  linearData[MHDTerm::BY] = _avZ[5];
  linearData[MHDTerm::BZ] = _avZ[6];
  linearData[MHDTerm::B] = sqrt(linearData[MHDTerm::BX]*linearData[MHDTerm::BX] +
				linearData[MHDTerm::BY]*linearData[MHDTerm::BY] +
				linearData[MHDTerm::BZ]*linearData[MHDTerm::BZ]);

  linearData[MHDTerm::P] = (_model->getGamma() - 1.0)*
    (_avZ[7] - 0.5*(linearData[MHDTerm::RHO]*linearData[MHDTerm::V]*
		    linearData[MHDTerm::V] + linearData[MHDTerm::B]*
		    linearData[MHDTerm::B]));
  linearData[MHDTerm::A] = sqrt(_model->getGamma()*linearData[MHDTerm::P]/
				linearData[MHDTerm::RHO]);
  
}

//////////////////////////////////////////////////////////////////////////////

} // namespace MHD
} // namespace Physics
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
