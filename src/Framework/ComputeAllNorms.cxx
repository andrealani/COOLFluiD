// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/NotImplementedException.hh"
#include "MathTools/MathConsts.hh"

#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"

#include "Framework/Storage.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/State.hh"
#include "Framework/DataHandle.hh"
#include "Framework/MeshData.hh"
#include "Framework/Framework.hh"
#include "Framework/ComputeAllNorms.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::MathTools;

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ComputeAllNorms, ComputeNorm, FrameworkLib, 1>
aComputeAllNormsProvider("AllNorms");

//////////////////////////////////////////////////////////////////////////////

ComputeAllNorms::ComputeAllNorms(const std::string& name) :
ComputeNorm(name),
m_gr(*this),
socket_rhs("rhs"),
socket_states("states")
{
}

//////////////////////////////////////////////////////////////////////////////

ComputeAllNorms::~ComputeAllNorms()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeAllNorms::GR_Combine (const CFreal & S1,
        const CFreal & S2,
        CFreal& D)
{
  throw Common::NotImplementedException (FromHere(),"ComputeAllNorms is not working\n");
  D = S1 + S2;
}

//////////////////////////////////////////////////////////////////////////////

CFreal ComputeAllNorms::GR_GetLocalValue () const
 {
    throw Common::NotImplementedException (FromHere(),"ComputeAllNorms is not working\n");

    cf_assert(m_var < PhysicalModelStack::getActive()->getNbEq());

    DataHandle< CFreal> rhs = socket_rhs.getDataHandle();
    DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
    const CFuint nbStates = rhs.size()/nbEqs;

    cf_assert (states.size()== nbStates);

    CFreal value = 0.0;
    for (CFuint i = 0; i < nbStates; ++i) {
      if ((states[i])->isParUpdatable()) {
        const CFreal tmp = rhs(i, m_compute_var_id[m_var_itr], nbEqs);
        value += tmp*tmp;
      }
    }

    return value;
 }

//////////////////////////////////////////////////////////////////////////////

RealVector ComputeAllNorms::compute ()
{
  throw Common::NotImplementedException (FromHere(),"ComputeAllNorms is not working\n");
  
  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  
  for(m_var_itr = 0; m_var_itr < m_residuals.size(); m_var_itr++){
    const CFreal globalValue = m_gr.GetGlobalValue (nsp);
    if(globalValue > 0.) {
      m_residuals[m_var_itr] = log10(sqrt(globalValue));
    }
    else {
      m_residuals[m_var_itr] = -MathTools::MathConsts::CFrealMax();
      //   m_residuals[m_var_itr] = log10(sqrt(MathTools::MathConsts::CFrealMin()));
    }
  }
  
  return m_residuals;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
ComputeAllNorms::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result = ComputeNorm::needsSockets();

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
