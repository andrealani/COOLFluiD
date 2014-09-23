// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MathTools/MathConsts.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/Framework.hh"
#include "Framework/SERComputeCFL.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<SERComputeCFL,
               ComputeCFL,
               FrameworkLib,
               1>
serComputeCFLProvider("SER");

//////////////////////////////////////////////////////////////////////////////

void SERComputeCFL::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal,Config::DynamicOption<> >("coeffCFL","Coefficient for the CFL calculation.");
   options.addConfigOption< CFreal,Config::DynamicOption<> >("power","power");
   options.addConfigOption< CFreal,Config::DynamicOption<> >("maxCFL","Max CFL.");
   options.addConfigOption< CFreal,Config::DynamicOption<> >("LimitCFL","Minimum value of the CFL");
   options.addConfigOption< bool,Config::DynamicOption<> >("Tol","Minimum value tolerence");
}

//////////////////////////////////////////////////////////////////////////////

SERComputeCFL::SERComputeCFL(const std::string& name) :
  ComputeCFL(name),
  m_Min(10.0),
  m_initialCFL(0.1),
  _Res0(10.0)
 {
  addConfigOptionsTo(this);

  
  m_coeffCFL = 1.001;
  setParameter("coeffCFL",&m_coeffCFL);

  m_maxCFL = 100;
  setParameter("maxCFL",&m_maxCFL);

  m_powerExp = 0.4;
  setParameter("power",&m_powerExp);
  
  m_LimitCFL = 4.;
  setParameter("LimitCFL",&m_LimitCFL);
  
  m_Tol = false;
  setParameter("Tol",&m_Tol);
 }

//////////////////////////////////////////////////////////////////////////////

SERComputeCFL::~SERComputeCFL()
{
}

//////////////////////////////////////////////////////////////////////////////

void SERComputeCFL::configure ( Config::ConfigArgs& args )
{
  ComputeCFL::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void SERComputeCFL::operator() (const ConvergenceStatus& m_cstatus)
{
  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();

   const CFuint iter = subSysStatus->getNbIter();
  if ( iter <= 5 )
   { 
        _Res0 = subSysStatus->getResidual();
        m_Min        = _Res0; 
        m_initialCFL = _cfl->getCFLValue(); 
   }
  else
    {
      const CFreal new_res = subSysStatus->getResidual();
      bool condstop = ((new_res<15)||(std::abs(new_res)<1e40));
      assert(condstop); 
      const CFreal Limiter = std::abs(new_res -_Res0);
            CFreal m_Max = (new_res >= -1.) ? m_Min + 1.0 : m_Min + 0.2 ;
      const CFreal ratio = _Res0/new_res;
          if ((Limiter >= 1.5) || ((new_res-m_Min)>=1.)) {
               m_Min = new_res;
           }      
           if (((new_res<= m_Min) || (new_res >= m_Max) || ((m_Tol)&&(new_res< _Res0))) && (ratio > 0) && (Limiter < 0.2))
               {      
              if (std::abs(new_res) > MathTools::MathConsts::CFrealEps())
              {
               const CFreal currentcfl = _cfl->getCFLValue();
               const CFreal limit_ratio = std::min<CFreal>(ratio,1.01);
               const CFreal invlimit_ratio = std::min<CFreal>(1./ratio,1.01);
                 CFreal cflvalue = (new_res <=0) ? m_coeffCFL*currentcfl*pow(invlimit_ratio, m_powerExp) : m_coeffCFL*currentcfl*pow(limit_ratio, m_powerExp);
               const CFreal minCFL = std::max(m_initialCFL,cflvalue);  
               const CFreal newCFL = (new_res > -1.) ? std::min(minCFL,m_LimitCFL): std::min(cflvalue,m_maxCFL);  
               _cfl->setCFLValue(newCFL);
                if ((new_res<=m_Min)){
               m_Min = new_res;
                 }     
              }
           }
              _Res0= new_res;
    }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
