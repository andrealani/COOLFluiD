// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/ComputeNorm.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/BadValueException.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void ComputeNorm::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("MonitoredVarID","ID of the variable whose residual will be monitored.");
   options.addConfigOption< std::vector<CFuint> >("ComputedVarID","IDs of the variables whose residual will be computed.");
   options.addConfigOption< bool >("NormalizedRes","Normalize the residual with the values");
   options.addConfigOption< std::vector<CFreal> >("RefVals","Values to normalize");
   options.addConfigOption< bool >("GlobalRes","Flag to use the global residual in the convergence method"); 
}

//////////////////////////////////////////////////////////////////////////////

ComputeNorm::ComputeNorm(const std::string & name) :
  NumericalStrategy(name),
  m_var_idx(0),
  m_var(0)
{
  addConfigOptionsTo(this);
  setParameter("MonitoredVarID",&m_var);

  m_compute_var_id = std::vector<CFuint>();
  setParameter("ComputedVarID",&m_compute_var_id);

  m_normalizedRes = false;
  setParameter("NormalizedRes",&m_normalizedRes);

  m_refVals = std::vector<CFreal>(); 
  setParameter("RefVals",&m_refVals);

  m_global_res = false;
  setParameter("GlobalRes",&m_global_res);
}


//////////////////////////////////////////////////////////////////////////////

ComputeNorm::~ComputeNorm()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeNorm::setup()
{
  CFAUTOTRACE;

  NumericalStrategy::setup();

  const CFuint nbEq = PhysicalModelStack::getActive()->getNbEq();

  if( !m_compute_var_id.empty() )
  {
    bool found = false;
    for(CFuint iVar=0; iVar < m_compute_var_id.size(); ++iVar)
    {
      if (m_compute_var_id[iVar] >= nbEq)
        throw Common::BadValueException (FromHere(),"ID of variable to compute the norm must be less than number of equations");

      if (m_compute_var_id[iVar] == m_var)
      {
        m_var_idx = iVar;
        found = true;
      }
    }

    cf_assert(found == true);
  }
  else // fill the list of variables if user did not
  {
    m_compute_var_id.resize(nbEq);
    for (CFuint i = 0; i < nbEq; ++i)
      m_compute_var_id[i] = i;
    m_var_idx = m_var;
  }

  m_residuals.resize(m_compute_var_id.size());
  
  m_refVals.resize(nbEq, 1.);
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSink> >
ComputeNorm::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  return result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
