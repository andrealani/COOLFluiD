// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NewtonMethod/NewtonIteratorData.hh"
#include "NewtonMethod/NewtonMethod.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace NewtonMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<NewtonIteratorData>,
		      NewtonIteratorData, NewtonMethodModule>
nullNewtonIteratorComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void NewtonIteratorData::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >        ("Norm","L2 Norm of dU to reach to stop the newton loop.");
   options.addConfigOption< bool >          ("SaveSystemToFile","Save files of matrix rhs solution vectors at each Newton step");
   options.addConfigOption< bool >          ("PrintHistory","Print convergence history for each Newton Iterator step");
   options.addConfigOption< vector<CFuint> >("MaxSteps","Maximum steps to perform in the newton loop.");
}

//////////////////////////////////////////////////////////////////////////////

NewtonIteratorData::NewtonIteratorData(Common::SafePtr<Framework::Method> owner)
  : ConvergenceMethodData(owner),
    m_achieved(false),
    m_lss()
{
   addConfigOptionsTo(this);

   m_maxSteps = vector<CFuint>();
   setParameter("MaxSteps",&m_maxSteps);

   m_maxNorm = -10.;
   setParameter("Norm",&m_maxNorm);

   m_printHistory = false;
   setParameter("PrintHistory",&m_printHistory);

  m_saveSystemToFile = false;
  setParameter("SaveSystemToFile",&m_saveSystemToFile);
}

//////////////////////////////////////////////////////////////////////////////

NewtonIteratorData::~NewtonIteratorData()
{
}

//////////////////////////////////////////////////////////////////////////////

void NewtonIteratorData::configure ( Config::ConfigArgs& args )
{
  ConvergenceMethodData::configure(args);

  // if the maximum number of steps has not been specified, just resize
  // the corresponding vector and set it to 1
  if (m_maxSteps.size() == 0) {
    m_maxSteps.resize(1);
    m_maxSteps[0] = 1;
  }
  cf_assert(m_maxSteps.size() > 0);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NewtonMethod

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
