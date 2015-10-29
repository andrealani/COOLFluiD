// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/ExprComputeDT.hh"
#include "Common/ParserException.hh"
#include "Framework/SubSystemStatus.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ExprComputeDT,
         ComputeDT,FrameworkLib,
         1>
exprComputeDTProvider("FunctionDT");

//////////////////////////////////////////////////////////////////////////////

void ExprComputeDT::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("Def","Definition of the Function.");
}

//////////////////////////////////////////////////////////////////////////////

ExprComputeDT::ExprComputeDT(const std::string& name) :
  ComputeDT(name),
  _functionParser(),
  _firstResidual(0.0),
  _lastResidual(0.0),
  _maxResidual(0.0),
  _eval(1.0,6),
  _vars("i,r,ri,rl,rmax,t")
{
   addConfigOptionsTo(this);
  _function = "1.0";
   setParameter("Def",&_function);
}

//////////////////////////////////////////////////////////////////////////////

ExprComputeDT::~ExprComputeDT()
{
}

//////////////////////////////////////////////////////////////////////////////

void ExprComputeDT::operator() ()
{
  const CFuint nbIter = SubSystemStatusStack::getActive()->getNbIter();
  // store the first residual
  if(nbIter == 1) {
    _firstResidual = SubSystemStatusStack::getActive()->getResidual();
  }

  _eval[0] = static_cast<CFreal>(nbIter);
  _eval[1] = SubSystemStatusStack::getActive()->getResidual();
  _eval[2] = _firstResidual;
  _eval[3] = _lastResidual;
  _eval[4] = _maxResidual;
  _eval[5] = SubSystemStatusStack::getActive()->getCurrentTimeDim();

  const CFreal currResidual = SubSystemStatusStack::getActive()->getResidual();

  _lastResidual = currResidual;
  if(_maxResidual < currResidual) {
    _maxResidual = currResidual;
  }

  if (nbIter > 1) {
    const CFreal DTnumber = _functionParser.Eval(&_eval[0]);// this is a dimensional time step!
    cf_assert(std::abs(DTnumber) > MathTools::MathConsts::CFrealEps());
    SubSystemStatusStack::getActive()->setDTDim(DTnumber);
    if (SubSystemStatusStack::getActive()->getTimeStepLayers() > 1){
      for (CFuint j=0; j < SubSystemStatusStack::getActive()->getTimeStepLayers() ; ++j){
        SubSystemStatusStack::getActive()->setInnerDT(j,SubSystemStatusStack::getActive()->getInnerDTRatio(j)*SubSystemStatusStack::getActive()->getDT());
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void ExprComputeDT::configure ( Config::ConfigArgs& args )
{
  ComputeDT::configure(args);

  try {
    setFunction();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }
}

//////////////////////////////////////////////////////////////////////////////

void ExprComputeDT::setFunction()
{
  // some sanity checks
  std::vector<std::string> _functionDef = Common::StringOps::getWords(_function);
  cf_assert(_functionDef.size() == 1);

  _functionParser.Parse(_function, _vars);

  if (_functionParser.ErrorMsg() != 0) {
    std::string msg("ParseError in DT::setFuntion(): ");
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
