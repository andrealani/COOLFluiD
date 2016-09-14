#include "ComputeDiffusiveFlux.hh"
#include "Framework/NamespaceSwitcher.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

ComputeDiffusiveFlux::ComputeDiffusiveFlux(const std::string& name) :
  ComputeFlux<CellCenterFVMData>(name),
  socket_normals("normals"), 
  socket_faceAreas("faceAreas"),
  socket_updateCoeff("updateCoeff"),
  socket_isOutward("isOutward"),
  socket_nstates("nstates"),
  _lFluxJacobian(),
  _rFluxJacobian()
{
}
      
//////////////////////////////////////////////////////////////////////////////

ComputeDiffusiveFlux::~ComputeDiffusiveFlux()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeDiffusiveFlux::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > ComputeDiffusiveFlux::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result  = ComputeFlux<CellCenterFVMData>::needsSockets();
  
  result.push_back(&socket_normals);
  result.push_back(&socket_faceAreas);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_isOutward);
  result.push_back(&socket_nstates);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void ComputeDiffusiveFlux::setup()
{
  ComputeFlux<CellCenterFVMData>::setup();
  
  _lFluxJacobian.resize(PhysicalModelStack::getActive()->getNbEq(),
			PhysicalModelStack::getActive()->getNbEq());
  
  _rFluxJacobian.resize(PhysicalModelStack::getActive()->getNbEq(),
			PhysicalModelStack::getActive()->getNbEq());
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
