// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"
#include "MathTools/MathConsts.hh"

#include "Framework/Storage.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/State.hh"
#include "Framework/DataHandle.hh"
#include "Framework/MeshData.hh"
#include "Framework/Framework.hh"
#include "Framework/ComputeL2Norm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::MathTools;

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ComputeL2Norm, ComputeNorm, FrameworkLib, 1>
computeL2NormProvider("L2");

//////////////////////////////////////////////////////////////////////////////

void ComputeL2Norm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("VectorName","Name of the vector of which to take the norm.");
  options.addConfigOption< CFreal >("Tolerance","If residual < tolerance residual is set to 0.");
}

//////////////////////////////////////////////////////////////////////////////

ComputeL2Norm::ComputeL2Norm(const std::string& name) :
ComputeNorm(name),
m_gr(*this),
sockets_norm(),
socket_states("states"),
m_vecnorm_name()
{
  addConfigOptionsTo(this);
  m_vecnorm_name = "rhs";
  setParameter("VectorName",&m_vecnorm_name);
  
  m_tolerance = 0.;
  setParameter("Tolerance",&m_tolerance);
}

//////////////////////////////////////////////////////////////////////////////

 ComputeL2Norm::~ComputeL2Norm()
 {
   CFLog(VERBOSE, "ComputeL2Norm::~ComputeL2Norm()\n");
 }

//////////////////////////////////////////////////////////////////////////////

void ComputeL2Norm::configure ( Config::ConfigArgs& args )
{
  ComputeNorm::configure(args);
  sockets_norm.createSocketSink<CFreal>(m_vecnorm_name);
}

//////////////////////////////////////////////////////////////////////////////

void ComputeL2Norm::GR_Combine (const CFreal & S1, const CFreal & S2, CFreal& D)
{
  D = S1 + S2;
}

//////////////////////////////////////////////////////////////////////////////

CFreal ComputeL2Norm::GR_GetLocalValue () const
{
  cf_assert(m_var < PhysicalModelStack::getActive()->getNbEq());
    
  DataHandle< CFreal > vecnorm = sockets_norm.getSocketSink<CFreal>(m_vecnorm_name)->getDataHandle();
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  
  CFLog(VERBOSE, "ComputeL2Norm::GR_GetLocalValue() => m_vecnorm_name = " << m_vecnorm_name << "\n");
  
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbStates = vecnorm.size()/nbEqs;
  
  CFLog(VERBOSE, "ComputeL2Norm::GR_GetLocalValue() => nbEqs = " << nbEqs << "\n");
  CFLog(VERBOSE, "ComputeL2Norm::GR_GetLocalValue() => vecnorm.size() = " << vecnorm.size() << "\n");
  
  if (states.size() != nbStates) {
    CFLog(VERBOSE, "ComputeL2Norm::GR_GetLocalValue() => " << states.size()  << " != " <<  nbStates << "\n");
    cf_assert (states.size()== nbStates);
  }
  
  const CFuint iVar = m_compute_var_id[m_var_itr];
  
  CFreal value = 0.0;
  for (CFuint i = 0; i < nbStates; ++i) {
    if ((states[i])->isParUpdatable()) {
      
      if (m_tolerance > 0.) {
	if (std::abs(vecnorm(i, iVar, nbEqs) < m_tolerance)) {
	  vecnorm(i, iVar, nbEqs) = 0.;
	}
      }
      
      const CFreal tmp = vecnorm(i, iVar, nbEqs);
      value += tmp*tmp; 
    }
  } 
  
  return value; 
}

//////////////////////////////////////////////////////////////////////////////

RealVector ComputeL2Norm::compute ()
{
  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();

  for(m_var_itr = 0; m_var_itr < m_residuals.size(); m_var_itr++)
  {
    CFreal globalValue = m_gr.GetGlobalValue (nsp);
    if(m_normalizedRes) { globalValue = globalValue/m_refVals[m_var_itr]; }
    if(globalValue > 0.)
    {
      m_residuals[m_var_itr] = log10(sqrt(globalValue));
    }
    else
    {
      m_residuals[m_var_itr] = -MathTools::MathConsts::CFrealMax();
    }
  }

  return m_residuals;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
ComputeL2Norm::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = ComputeNorm::needsSockets();
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > sockets = sockets_norm.getAllSinkSockets();

  result.insert( result.end(), sockets.begin(), sockets.end() );

  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
