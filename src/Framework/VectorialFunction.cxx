// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/CFLog.hh"

#include "VectorialFunction.hh"
#include "MathTools/FunctionParser.hh"
#include <math.h>
//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

VectorialFunction::VectorialFunction()
  : m_isParsed(false),
    m_vars(""),
    m_nbVars(0),
    m_functions(0),
    m_parsers(),
    m_result()
{
}

//////////////////////////////////////////////////////////////////////////////

VectorialFunction::~VectorialFunction()
{
  clear();
}

//////////////////////////////////////////////////////////////////////////////

void VectorialFunction::clear()
{
  m_isParsed = false;
  m_result.resize(0);
  for(CFuint i = 0; i < m_parsers.size(); i++) {
      deletePtr(m_parsers[i]);
  }
  vector<FunctionParser*>().swap(m_parsers);
}

//////////////////////////////////////////////////////////////////////////////

void VectorialFunction::setVariables(const vector<std::string>& vars)
{
  m_nbVars = vars.size();
  m_vars = "";
  if(vars.size() > 0) {
    m_vars = vars[0];
    for(CFuint i = 1; i < m_nbVars; i++) {
      m_vars += ",";
      m_vars += vars[i];
    }
  }
  
  CFLog(DEBUG_MAX, "VectorialFunction::setVariables() => m_nbVars = "
	<< m_nbVars << ", vars = " << m_vars << "\n");
  
  m_isParsed = false;
}
    
//////////////////////////////////////////////////////////////////////////////

void VectorialFunction::setFunctions(const vector<std::string>& functions)
{
  m_functions = functions;
  m_isParsed = false;
}

//////////////////////////////////////////////////////////////////////////////

void VectorialFunction::parse()
{
  CFAUTOTRACE;

  clear();

  for(CFuint i = 0; i < m_functions.size(); i++) {

    FunctionParser* ptr = new FunctionParser();
    
    m_parsers.push_back(ptr);
    CFLogDebugMax("Parsing Function: " << m_functions[i] << " Vars: " << m_vars << "\n");
    ptr->Parse(m_functions[i],m_vars);
    
    if (ptr->ErrorMsg() != 0) {
      std::string msg("ParseError in VectorialFunction::parse(): ");
      msg += std::string(ptr->ErrorMsg());
      msg += " Function: " + m_functions[i];
      msg += " Vars: "     + m_vars;
      throw Common::ParserException (FromHere(),msg);
    }
  }
  
  m_result.resize(m_functions.size());
  m_isParsed = true;
}

//////////////////////////////////////////////////////////////////////////////

void VectorialFunction::evaluate(const RealVector& varValue,
				 RealVector& value) const
{
  cf_assert(m_isParsed);
  if (varValue.size() != m_nbVars) { 
    CFLog(ERROR, "varValue.size()=[" << varValue.size() << "] != m_nbVars=[" << m_nbVars << "]\n");
    cf_assert(varValue.size() == m_nbVars);
  }

  // evaluate and store the functions line by line in the vector
  for(CFuint i = 0; i < m_parsers.size(); i++) {
    value[i] = m_parsers[i]->Eval(&const_cast<RealVector&>(varValue)[0]);
  }
}

//////////////////////////////////////////////////////////////////////////////

void VectorialFunction::evaluate(CFuint iVar, const RealVector& varValue, 
				 CFreal& value)
{
  cf_assert(m_isParsed);
  cf_assert(iVar < m_parsers.size());

  // evaluate and store the functions line by line in the vector
  value = m_parsers[iVar]->Eval(&const_cast<RealVector&>(varValue)[0]);
}

//////////////////////////////////////////////////////////////////////////////

RealVector& VectorialFunction::operator()(const RealVector& varValue)
{
  cf_assert(m_isParsed);
  cf_assert(varValue.size() == m_nbVars);

  // evaluate and store the functions line by line in the vector
  for(CFuint i = 0; i < m_parsers.size(); i++) {
    m_result[i] = m_parsers[i]->Eval(&const_cast<RealVector&>(varValue)[0]);
  }
  
  return m_result;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

