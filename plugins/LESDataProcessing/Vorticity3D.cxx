
#include "Vorticity3D.hh"
#include "Framework/DataHandle.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "LESDataProcessing.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
	
		namespace LESDataProcessing {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider<Vorticity3D, 
                       LESProcessingData, 
                       TurbulenceFunction, 
                       LESDataProcessingModule> 
vorticity3DProvider("Vorticity3D");

//////////////////////////////////////////////////////////////////////////////

void Vorticity3D::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

Vorticity3D::Vorticity3D(const std::string& name) :
  TurbulenceFunction(name)
{
  addConfigOptionsTo(this);
  // Setting default configurations here.  
}

//////////////////////////////////////////////////////////////////////////////


Vorticity3D::~Vorticity3D()
{
}

//////////////////////////////////////////////////////////////////////////////

void Vorticity3D::configure ( Config::ConfigArgs& args )
{
  TurbulenceFunction::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void Vorticity3D::setup()
{
  CFAUTOTRACE; 
  TurbulenceFunction::setup();
  // Set the size of the socket
  const CFuint dim = Framework::PhysicalModelStack::getActive()->getDim();
  m_socket.resize(getNbStates()*dim);  
}

//////////////////////////////////////////////////////////////////////////////

void Vorticity3D::compute(const RealVector& state, std::vector<RealVector>& grad, const CFuint& iCell)
{
  
  // assume gradients of primitive variables
  m_socket(iCell,XX,3) = grad[W][YY] - grad[V][ZZ];      //    dw/dy - dv/dz
  m_socket(iCell,YY,3) = grad[U][ZZ] - grad[W][XX];      //    du/dz - dw/dx
  m_socket(iCell,ZZ,3) = grad[V][XX] - grad[U][YY];      //    dv/dx - du/dy
  
}
//////////////////////////////////////////////////////////////////////////////
	
	  } // end of namespace LESDataProcessing
		
  } // namespace Numerics

} // namespace COOLFluiD
