
#include "VelocityGradients3D.hh"
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

Framework::MethodStrategyProvider<VelocityGradients3D, 
                       LESProcessingData, 
                       TurbulenceFunction, 
                       LESDataProcessingModule> 
velocityGradients3DProvider("VelocityGradients3D");

//////////////////////////////////////////////////////////////////////////////

void VelocityGradients3D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("VelocityComponents","Velocity components to store: (e.g. \"U V W\" for all)");
}

//////////////////////////////////////////////////////////////////////////////

VelocityGradients3D::VelocityGradients3D(const std::string& name) :
  TurbulenceFunction(name),
  m_velocityComponents()
{
  addConfigOptionsTo(this);
  // Setting default configurations here.  
  m_velocityComponents = std::vector<std::string>();
  setParameter("VelocityComponents",&m_velocityComponents);
}

//////////////////////////////////////////////////////////////////////////////


VelocityGradients3D::~VelocityGradients3D()
{
}

//////////////////////////////////////////////////////////////////////////////

void VelocityGradients3D::configure ( Config::ConfigArgs& args )
{
  TurbulenceFunction::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void VelocityGradients3D::setup()
{
  CFAUTOTRACE; 
  TurbulenceFunction::setup();
  
  if(m_velocityComponents.size()==0) 
    m_velocityComponents.push_back("U");
  
  // Set the size of the socket
  m_dim = Framework::PhysicalModelStack::getActive()->getDim();
  m_nbVelocityComponents = m_velocityComponents.size();
  m_socket.resize(getNbStates()*m_dim*m_nbVelocityComponents);  
  
  m_comp.resize(m_nbVelocityComponents);
  m_stride = m_dim*m_comp.size();

  for(CFuint i=0; i<m_nbVelocityComponents; ++i) {
    if (m_velocityComponents[i] == "U") {
      m_comp[i] = U;
    }
    if (m_velocityComponents[i] == "V") {
      m_comp[i] = V;
    }
    if (m_velocityComponents[i] == "W") {
      m_comp[i] = W;
    }
  }
  
}

//////////////////////////////////////////////////////////////////////////////

void VelocityGradients3D::compute(const RealVector& state, std::vector<RealVector>& grad, const CFuint& iCell)
{
  
  // assume gradients of primitive variables
  for(CFuint i=0; i<m_comp.size(); ++i) {
    m_socket(iCell, XX +i*m_dim,m_stride) = grad[m_comp[i]][XX];
    m_socket(iCell, YY +i*m_dim,m_stride) = grad[m_comp[i]][YY];
    m_socket(iCell, ZZ +i*m_dim,m_stride) = grad[m_comp[i]][ZZ];
  }
}
//////////////////////////////////////////////////////////////////////////////
	
	  } // end of namespace LESDataProcessing
		
  } // namespace Numerics

} // namespace COOLFluiD
