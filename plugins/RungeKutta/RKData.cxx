// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RungeKutta/RungeKutta.hh"

#include "RKData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RungeKutta {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<RKData>, RKData, RungeKuttaModule> nullRKComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void RKData::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< bool >("TimeAccurate","True if time accurate");
   options.addConfigOption< CFuint >("Order","Order of the R-K method");
   options.addConfigOption< std::vector<CFreal> >("Alpha","Alpha Coeficients needed for the R-K steps.");
   options.addConfigOption< std::vector<CFreal> >("Beta","Beta Coeficients needed for the R-K steps.");
   options.addConfigOption< std::vector<CFreal> >("Gamma","Gamma (time) Coeficients needed for the R-K steps.");
}

//////////////////////////////////////////////////////////////////////////////

RKData::RKData(Common::SafePtr<Framework::Method> owner)
  : ConvergenceMethodData(owner),
    m_isFirst(true),
    m_isLast(false),
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

RKData::~RKData()
{
}

//////////////////////////////////////////////////////////////////////////////

void RKData::configure ( Config::ConfigArgs& args )
{
  ConvergenceMethodData::configure(args);

  if(m_order==0)
  {
    CFout << "Order of RungeKutta Method not Defined " << "\n";
    CFout << "By Default: Using 4th order Method " << "\n";
    m_order = 4;
  }

  if((m_alpha.size() != m_order) || (m_beta.size() != m_order)) {
    CFout << "RK Coef Values Not Valid" << "\n";

    // Set defaults Values
    if (m_order == 2)
    {
      CFout << "Using 2nd order R-K with defaults coeficients" << "\n";
      m_alpha.resize(m_order);
      m_beta .resize(m_order);
      m_gamma.resize(m_order);

      // Coeficient from Metodos Numericos, Heitor Pina, pp516 (Euler-Cauchy Method)
      m_alpha[0]=.0;
      m_alpha[1]=.5;

      m_beta[0]=0.;
      m_beta[1]=1.;

      m_gamma[0] = 0.0;
      m_gamma[1] = 0.5;
    }
    else {
      if (m_order == 3) {
        CFout << "Using 3rd Order RK with defaults coeficients" << "\n";
        m_alpha.resize(m_order);
        m_beta .resize(m_order);
        m_gamma.resize(m_order);

        // Coeficient from Metodos Numericos, Heitor Pina, pp516 (Heun Method)
        m_alpha[0] = 0.;
        m_alpha[1] = 1./6.;
        m_alpha[2] = 2./6.;

        m_beta[0] = 0.25;
        m_beta[1] = 0.;
        m_beta[2] = 0.75;

        m_gamma[0] = 0.0;
        m_gamma[1] = 0.5;
        m_gamma[2] = 1.0;
      }
      else
      {
        CFout << "Using 4th Order RK with defaults coeficients" << "\n";
        m_order = 4;
        m_alpha.resize(m_order);
        m_beta .resize(m_order);
        m_gamma.resize(m_order);

        // Coeficient from Numerical Computation of Internal and External, C. Hirsch, pp446
        m_alpha[0] = 0.;
        m_alpha[1] = .5;
        m_alpha[2] = .5;
        m_alpha[3] = 1.;

        m_beta[0] = 1./6.;
        m_beta[1] = 1./3.;
        m_beta[2] = 1./3.;
        m_beta[3] = 1./6.;


        m_gamma[0] = 0.0;
        m_gamma[1] = 0.5;
        m_gamma[2] = 0.5;
        m_gamma[3] = 1.0;
      }
    }
  }
  //for (CFuint i = 0; i < m_order; i++){
  //CFout << "Alpha " << i << " = " << m_alpha[i] << "\n";
  //CFout << "Beta " << i << " = " << m_beta[i] << "\n";
  //}
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RungeKutta

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

