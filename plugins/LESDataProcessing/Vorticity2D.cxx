
#include "Vorticity2D.hh"
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

Framework::MethodStrategyProvider<Vorticity2D, 
                       LESProcessingData, 
                       TurbulenceFunction, 
                       LESDataProcessingModule> 
vorticity2DProvider("Vorticity2D");

//////////////////////////////////////////////////////////////////////////////

void Vorticity2D::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

Vorticity2D::Vorticity2D(const std::string& name) :
  TurbulenceFunction(name)
{
  addConfigOptionsTo(this);
  // Setting default configurations here.  
}

//////////////////////////////////////////////////////////////////////////////


Vorticity2D::~Vorticity2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void Vorticity2D::configure ( Config::ConfigArgs& args )
{
  TurbulenceFunction::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void Vorticity2D::setup()
{
  CFAUTOTRACE; 
  TurbulenceFunction::setup();
  // Set the size of the socket
  m_socket.resize(getNbStates());  
}

//////////////////////////////////////////////////////////////////////////////

void Vorticity2D::compute(const RealVector& state, std::vector<RealVector>& grad, const CFuint& iCell)
{
  
  // assume gradients of primitive variables
  m_socket[iCell] = grad[V][XX] - grad[U][YY];      //    dv/dx - du/dy
  
}
//////////////////////////////////////////////////////////////////////////////
	
	  } // end of namespace LESDataProcessing
		
  } // namespace Numerics

} // namespace COOLFluiD
