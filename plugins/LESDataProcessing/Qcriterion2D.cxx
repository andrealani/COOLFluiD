
#include "Qcriterion2D.hh"
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

Framework::MethodStrategyProvider<Qcriterion2D, 
                       LESProcessingData, 
                       TurbulenceFunction, 
                       LESDataProcessingModule> 
Qcriterion2DProvider("Qcriterion2D");

//////////////////////////////////////////////////////////////////////////////

void Qcriterion2D::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

Qcriterion2D::Qcriterion2D(const std::string& name) :
  TurbulenceFunction(name)
{
  addConfigOptionsTo(this);
  // Setting default configurations here.  
}

//////////////////////////////////////////////////////////////////////////////


Qcriterion2D::~Qcriterion2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void Qcriterion2D::configure ( Config::ConfigArgs& args )
{
  TurbulenceFunction::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void Qcriterion2D::setup()
{
  CFAUTOTRACE; 
  TurbulenceFunction::setup();

  // Set the size of the socket
  m_socket.resize(getNbStates());  
}

//////////////////////////////////////////////////////////////////////////////

void Qcriterion2D::compute(const RealVector& state, std::vector<RealVector>& gradients, const CFuint& iCell)
{  
  CFreal vorticity = getAbsoluteVorticity(gradients);
  CFreal strainRate = getAbsoluteStrainRate(gradients);
  m_socket[iCell] = 0.25*(vorticity*vorticity - strainRate*strainRate);  
}

//////////////////////////////////////////////////////////////////////////////

CFreal Qcriterion2D::getAbsoluteVorticity(const std::vector<RealVector>& grad) const {
  
  // assume gradients of primitive variables
  
  return grad[V][XX] - grad[U][YY];      //    dv/dx - du/dy

}

//////////////////////////////////////////////////////////////////////////////

CFreal Qcriterion2D::getAbsoluteStrainRate(const std::vector<RealVector>& grad) const {
    
    // Strain Rate Tensor
    RealMatrix S(2,2,0.0);
    S(XX,XX) = grad[U][XX];
    S(YY,YY) = grad[V][YY];
    S(XX,YY) = S(YY,XX) = 0.5*(grad[U][YY] + grad[V][XX]);
    
    // S[i,j]*S[i,j]
    CFreal SS = S(XX,XX)*S(XX,XX) + S(YY,YY)*S(YY,YY) + 2.0*(S(XX,YY)*S(XX,YY));
    // sqrt(2*S*S)
    return std::sqrt(2.0*SS);
}

//////////////////////////////////////////////////////////////////////////////
	
	  } // end of namespace LESDataProcessing
		
  } // namespace Numerics

} // namespace COOLFluiD
