// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/LSSData.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void LSSData::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint , Config::DynamicOption<> >("MaxIter","Maximum number of steps to be taken in the iterative solver.");
  options.addConfigOption< CFuint , Config::DynamicOption<> >("SaveRate","Save rate to write matrix and vector to file.");
  options.addConfigOption< CFuint , Config::DynamicOption<> >("PreconditionerRate","Rate at which recomputing the preconditioner.");
  options.addConfigOption< bool >("SaveSystemToFile","Save files of matrix rhs solution vectors at each solve");
  options.addConfigOption< bool >("Output","Flag indicating to output or not the solver convergence history.");
  options.addConfigOption< bool >("UseBlockPreconditioner","Enables to use block diagonal preconditioner (Dimension of each block is the number of local states on each CPU).");
  options.addConfigOption< bool >("UseGPU","Enables to run on GPU.");
  options.addConfigOption< bool >("UseNodeBased","Use node-based sparsity and assembly.");
}

//////////////////////////////////////////////////////////////////////////////

LSSData::LSSData(SafePtr<std::valarray<bool> > maskArray,
    CFuint& nbSysEquations,
    Common::SafePtr<Method> owner) : MethodData(owner),
    m_localToGlobal(),
    m_localToLocallyUpdateble(),
    m_maskArray(maskArray),
    m_nbSysEquations(nbSysEquations)
{
  addConfigOptionsTo(this);
  cf_assert(maskArray.isNotNull());

  m_maxIter = 50;
  setParameter("MaxIter",&m_maxIter);
  
  m_saveRate = 0;
  setParameter("SaveRate",&m_saveRate);

  m_preconditionerRate = 1; 
  setParameter("PreconditionerRate", &m_preconditionerRate);
  
  m_isOutput = false;
  setParameter("Output",&m_isOutput);
  
  m_saveSystemToFile = false;
  setParameter("SaveSystemToFile",&m_saveSystemToFile);
  
  _useBlockPreconditionerMatrix = false;
  setParameter("UseBlockPreconditioner", &_useBlockPreconditionerMatrix);
  
  m_useGPU = false;
  setParameter("UseGPU", &m_useGPU);
  
  m_useNodeBased = false;
  setParameter("UseNodeBased", &m_useNodeBased);
}

//////////////////////////////////////////////////////////////////////////////

LSSData::~LSSData()
{
}

//////////////////////////////////////////////////////////////////////////////

void LSSData::configure ( Config::ConfigArgs& args )
{
  MethodData::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void LSSData::setup()
{
  MethodData::setup();
}

//////////////////////////////////////////////////////////////////////////////

void LSSData::unsetup()
{
  MethodData::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
