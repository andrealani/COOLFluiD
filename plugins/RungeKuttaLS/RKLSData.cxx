// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RungeKuttaLS/RungeKuttaLS.hh"

#include "RKLSData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RungeKuttaLS {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<RKLSData>, RKLSData, RungeKuttaLSModule> nullRKLSComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void RKLSData::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< bool >("TimeAccurate","True if time accurate");
   options.addConfigOption< CFuint >("Order","Order of the R-K method");
   options.addConfigOption< std::vector<CFreal> >("Alpha","Alpha Coeficients needed for the R-K steps.");
   options.addConfigOption< std::vector<CFreal> >("Beta","Beta Coeficients needed for the R-K steps.");
   options.addConfigOption< std::vector<CFreal> >("Gamma","Gamma (time) Coeficients needed for the R-K steps.");
}

//////////////////////////////////////////////////////////////////////////////

RKLSData::RKLSData(Common::SafePtr<Framework::Method> owner)
  : ConvergenceMethodData(owner),
    m_alpha(0),
    m_beta(0),
    m_gamma(0),
    m_timeIterN()
{
  addConfigOptionsTo(this);

  m_order = 4;
  setParameter("Order",&m_order);

  setParameter("Alpha",&m_alpha);

  setParameter("Beta" ,&m_beta );

  setParameter("Gamma",&m_gamma);

  m_isTimeAccurate = false;
  setParameter("TimeAccurate",&m_isTimeAccurate);
}

//////////////////////////////////////////////////////////////////////////////

RKLSData::~RKLSData()
{
}

//////////////////////////////////////////////////////////////////////////////

void RKLSData::configure ( Config::ConfigArgs& args )
{
  ConvergenceMethodData::configure(args);

  if(m_order == 0)
  {
    CFout << "Order of RungeKuttaLS Method not Defined " << "\n";
    CFout << "By Default: Using 4th order Method " << "\n";
    m_order = 4;
  }

  if((m_alpha.size() != m_order) || (m_beta.size() != m_order) || (m_gamma.size() != m_order))
  {
    CFout << "R-K Coef Values Not Valid" << "\n";

    // Set defaults Values
    if (m_order == 2)
    {
      CFout << "Using 2nd order R-K with defaults coeficients" << "\n";
      m_alpha.resize(m_order);
      m_beta .resize(m_order);
      m_gamma.resize(m_order);

      m_alpha[0] = 0.0;
      m_alpha[1] = 0.0;

      m_beta[0] = 0.5;
      m_beta[1] = 1.0;

      m_gamma[0] = 0.0;
      m_gamma[1] = 0.5;
    }
    else if (m_order == 3)
    {
        CFout << "Using 3rd order R-K with defaults coeficients" << "\n";
        m_alpha.resize(m_order);
        m_beta .resize(m_order);
        m_gamma.resize(m_order);

        // 3rd order TVD R-K scheme
        m_alpha[0] = 0.0;
        m_alpha[1] = 1.0/4.0;
        m_alpha[2] = 2.0/3.0;

        m_beta[0] = 1.0;
        m_beta[1] = 1.0/4.0;
        m_beta[2] = 2.0/3.0;

        m_gamma[0] = 0.0;
        m_gamma[1] = 0.5;
        m_gamma[2] = 1.0;
    }
    else
    {
      CFout << "Using 4th order R-K with defaults coeficients" << "\n";
      m_order = 4;
      m_alpha.resize(m_order);
      m_beta .resize(m_order);
      m_gamma.resize(m_order);

      // 4th order scheme
      m_alpha[0] = 0.0;
      m_alpha[1] = 0.0;
      m_alpha[2] = 0.0;
      m_alpha[3] = 0.0;

      m_beta [0] = 1.0/4.0;
      m_beta [1] = 1.0/3.0;
      m_beta [2] = 1.0/2.0;
      m_beta [3] = 1.0;

      m_gamma[0] = 0.0;
      m_gamma[1] = 0.5;
      m_gamma[2] = 0.5;
      m_gamma[3] = 1.0;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RungeKuttaLS

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

