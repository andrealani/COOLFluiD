// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/ErrorEstimatorMethod.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void ErrorEstimatorMethod::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("EstimateRate","Rate to estimate the error in the solution.");
}

//////////////////////////////////////////////////////////////////////////////

void ErrorEstimatorMethod::build_dynamic_functions()
{
  add_dynamic_function("estimate",&ErrorEstimatorMethod::estimate);
}

//////////////////////////////////////////////////////////////////////////////

void ErrorEstimatorMethod::registActionListeners()
{
  Method::registActionListeners();
  
  // const std::string ssname = SubSystemStatusStack::getCurrentName();   
  //   event_handler->addListener(event_handler->key(ssname, "CF_ON_ERRORESTIMATOR_ESTIMATE"),
  //                              this,&ErrorEstimatorMethod::estimate);
}

//////////////////////////////////////////////////////////////////////////////

ErrorEstimatorMethod::ErrorEstimatorMethod(const std::string& name) : Method(name)
{
  // define which functions might be called dynamic
  build_dynamic_functions();
  // regist which functions might be called by raising Events
  registActionListeners();
  // regist the configuration options
  addConfigOptionsTo(this);

  _estimateRate = 10;
  setParameter("EstimateRate",&_estimateRate);
}

//////////////////////////////////////////////////////////////////////////////

ErrorEstimatorMethod::~ErrorEstimatorMethod()
{
}

//////////////////////////////////////////////////////////////////////////////

void ErrorEstimatorMethod::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  Method::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void ErrorEstimatorMethod::estimate()
{
  CFAUTOTRACE;

  cf_assert(isConfigured());
  cf_assert(isSetup());
  //cf_assert(isSpaceMethodSet());

  pushNamespace();

  estimateImpl();

  popNamespace();
}

//////////////////////////////////////////////////////////////////////////////

void ErrorEstimatorMethod::setMethodImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

void ErrorEstimatorMethod::unsetMethodImpl()
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD
