// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ObjectProvider.hh"
#include "Framework/ExprComputeCFL.hh"
#include "Common/ParserException.hh"
#include "Framework/CFL.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ExprComputeCFL, ComputeCFL,FrameworkLib, 1>
epxrComputeCFLProvider("Function");

//////////////////////////////////////////////////////////////////////////////

void ExprComputeCFL::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("Def","Definition of the Function.");
}

//////////////////////////////////////////////////////////////////////////////

ExprComputeCFL::ExprComputeCFL(const std::string& name) :
  ComputeCFL(name),
  _functionParser(),
  _firstResidual(0.0),
  _lastResidual(0.0),
  _maxResidual(0.0),
  _eval(1.0,7),
  _vars("i,r,ri,rl,rmax,cfl,si")
{
   addConfigOptionsTo(this);
  _function = "1.0";
   setParameter("Def",&_function);
}

//////////////////////////////////////////////////////////////////////////////

ExprComputeCFL::~ExprComputeCFL()
{
}

//////////////////////////////////////////////////////////////////////////////

void ExprComputeCFL::operator() (const ConvergenceStatus& cstatus)
{
  // store the first residual
  if(SubSystemStatusStack::getActive()->getNbIter() == 1) {
    _firstResidual = SubSystemStatusStack::getActive()->getResidual();
  }
  
  _eval[0] = cstatus.iter;
  _eval[1] = cstatus.res;
  _eval[2] = _firstResidual;
  _eval[3] = _lastResidual;
  _eval[4] = _maxResidual;
  _eval[5] = _cfl->getCFLValue();
  _eval[6] = cstatus.subiter;
  
  const CFreal currResidual = SubSystemStatusStack::getActive()->getResidual();
  
  _lastResidual = currResidual;
  if(_maxResidual < currResidual) {
    _maxResidual = currResidual;
  }
  
  //   if (cstatus.iter > 1) {
  const CFreal CFLnumber = _functionParser.Eval(&_eval[0]);
  CFLog(VERBOSE, "ExprComputeCFL::operator() => _eval = " << _eval 
	<< ", CFLnumber = " << CFLnumber << "\n");
  cf_assert(std::abs(CFLnumber) > MathTools::MathConsts::CFrealEps());
  _cfl->setCFLValue(CFLnumber);
  // }
}

//////////////////////////////////////////////////////////////////////////////

void ExprComputeCFL::configure ( Config::ConfigArgs& args )
{
  ComputeCFL::configure(args);

  try {
    setFunction();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }
}

//////////////////////////////////////////////////////////////////////////////

void ExprComputeCFL::setFunction()
{
  // some sanity checks
  std::vector<std::string> functionDef = Common::StringOps::getWords(_function);
  cf_assert(functionDef.size() == 1);

  _functionParser.Parse(_function, _vars);

  if (_functionParser.ErrorMsg() != 0) {
    std::string msg("ParseError in CFL::setFuntion(): ");
    msg += std::string(_functionParser.ErrorMsg());
    msg += " Function: " + _function;
    msg += " Vars: "     + _vars;
    throw Common::ParserException (FromHere(),msg);
  }
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
