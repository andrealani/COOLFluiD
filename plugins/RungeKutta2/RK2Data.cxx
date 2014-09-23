// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "RungeKutta2/RungeKutta2.hh"

#include "RK2Data.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace RungeKutta2 {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<RK2Data>, RK2Data, RungeKutta2Module> nullRK2ComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void RK2Data::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< bool >("TimeAccurate","True if time accurate");
}

//////////////////////////////////////////////////////////////////////////////

RK2Data::RK2Data(Common::SafePtr<Framework::Method> owner)
 : ConvergenceMethodData(owner)
{
  addConfigOptionsTo(this);

  m_isTimeAccurate = false;
  setParameter("TimeAccurate",&m_isTimeAccurate);
}

//////////////////////////////////////////////////////////////////////////////

RK2Data::~RK2Data()
{
}

//////////////////////////////////////////////////////////////////////////////

void RK2Data::configure ( Config::ConfigArgs& args )
{
  ConvergenceMethodData::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace RungeKutta2

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

