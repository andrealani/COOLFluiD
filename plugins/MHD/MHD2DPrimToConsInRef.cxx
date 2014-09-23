#include "MHD/MHD.hh"
#include "MHD2DPrimToConsInRef.hh"
#include "MHDPhysicalModel.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/State.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MHD2DPrimToConsInRef, VarSetMatrixTransformer, MHDModule, 1> mhd2DPrimToConsInRefProvider("MHD2DPrimToConsInRef");

//////////////////////////////////////////////////////////////////////////////

MHD2DPrimToConsInRef::MHD2DPrimToConsInRef(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
   _model(model->getConvectiveTerm().d_castTo<MHDTerm>())
{
}

//////////////////////////////////////////////////////////////////////////////

MHD2DPrimToConsInRef::~MHD2DPrimToConsInRef()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DPrimToConsInRef::setMatrixFromRef()
{
  const RealVector& linearData = _model->getPhysicalData();

  _transMatrix(0,0) = 1.0;

  _transMatrix(1,0) = linearData[MHDTerm::VX];
  _transMatrix(1,1) = linearData[MHDTerm::RHO];

  _transMatrix(2,0) = linearData[MHDTerm::VY];
  _transMatrix(2,2) = linearData[MHDTerm::RHO];

  _transMatrix(3,0) = linearData[MHDTerm::VZ];
  _transMatrix(3,3) = linearData[MHDTerm::RHO];

  _transMatrix(4,4) = 1.0;

  _transMatrix(5,5) = 1.0;

  _transMatrix(6,6) = 1.0;

  _transMatrix(7,0) = 0.5*(linearData[MHDTerm::V]*linearData[MHDTerm::V]);
  _transMatrix(7,1) = linearData[MHDTerm::RHO]*linearData[MHDTerm::VX];
  _transMatrix(7,2) = linearData[MHDTerm::RHO]*linearData[MHDTerm::VY];
  _transMatrix(7,3) = linearData[MHDTerm::RHO]*linearData[MHDTerm::VZ];
  _transMatrix(7,4) = linearData[MHDTerm::BX];
  _transMatrix(7,5) = linearData[MHDTerm::BY];
  _transMatrix(7,6) = linearData[MHDTerm::BZ];
  _transMatrix(7,7) = 1./(_model->getGamma() - 1.0);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
