
#include "Qcriterion3D.hh"
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

Framework::MethodStrategyProvider<Qcriterion3D, 
                       LESProcessingData, 
                       TurbulenceFunction, 
                       LESDataProcessingModule> 
Qcriterion3DProvider("Qcriterion3D");

//////////////////////////////////////////////////////////////////////////////

void Qcriterion3D::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

Qcriterion3D::Qcriterion3D(const std::string& name) :
  TurbulenceFunction(name)
{
  addConfigOptionsTo(this);
  // Setting default configurations here.  
}

//////////////////////////////////////////////////////////////////////////////


Qcriterion3D::~Qcriterion3D()
{
}

//////////////////////////////////////////////////////////////////////////////

void Qcriterion3D::configure ( Config::ConfigArgs& args )
{
  TurbulenceFunction::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void Qcriterion3D::setup()
{
  CFAUTOTRACE; 
  TurbulenceFunction::setup();

  // Set the size of the socket
  m_socket.resize(getNbStates());  
}

//////////////////////////////////////////////////////////////////////////////

void Qcriterion3D::compute(const RealVector& state, std::vector<RealVector>& gradients, const CFuint& iCell)
{  
  CFreal vorticity = getAbsoluteVorticity(gradients);
  CFreal strainRate = getAbsoluteStrainRate(gradients);
  m_socket[iCell] = 0.25*(vorticity*vorticity - strainRate*strainRate);  
}

//////////////////////////////////////////////////////////////////////////////

CFreal Qcriterion3D::getAbsoluteVorticity(const std::vector<RealVector>& grad) const {
  
  // assume gradients of primitive variables

  RealVector vorticity(3);
  vorticity[XX] = grad[W][YY] - grad[V][ZZ];      //    dw/dy - dv/dz
  vorticity[YY] = grad[U][ZZ] - grad[W][XX];      //    du/dz - dw/dx
  vorticity[ZZ] = grad[V][XX] - grad[U][YY];      //    dv/dx - du/dy
  
  return vorticity.norm2();

}

//////////////////////////////////////////////////////////////////////////////

CFreal Qcriterion3D::getAbsoluteStrainRate(const std::vector<RealVector>& grad) const {
    
    // Strain Rate Tensor
    RealMatrix S(3,3,0.0);
    S(XX,XX) = grad[U][XX];
    S(YY,YY) = grad[V][YY];
    S(ZZ,ZZ) = grad[W][ZZ];  
    S(XX,YY) = S(YY,XX) = 0.5*(grad[U][YY] + grad[V][XX]);
    S(XX,ZZ) = S(ZZ,XX) = 0.5*(grad[U][ZZ] + grad[W][XX]);
    S(YY,ZZ) = S(ZZ,YY) = 0.5*(grad[V][ZZ] + grad[W][YY]);
    
    // S[i,j]*S[i,j]
    CFreal SS = S(XX,XX)*S(XX,XX) + S(YY,YY)*S(YY,YY) + S(ZZ,ZZ)*S(ZZ,ZZ)
                + 2.0*(S(XX,YY)*S(XX,YY) + S(YY,ZZ)*S(YY,ZZ) + S(XX,ZZ)*S(XX,ZZ));
    // sqrt(2*S*S)
    return std::sqrt(2.0*SS);
}

//////////////////////////////////////////////////////////////////////////////
	
	  } // end of namespace LESDataProcessing
		
  } // namespace Numerics

} // namespace COOLFluiD
