// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/Framework.hh"
#include "Framework/LookupInterpolator.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"
#include "Common/NotImplementedException.hh"
#include "Common/PE.hh"
#include "Common/SwapEmpty.hh"
#include "Environment/DirPaths.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/ObjectProvider.hh"
#include "Environment/SingleBehaviorFactory.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

ObjectProvider<LookupInterpolator, StateInterpolator, FrameworkLib, 1>
lookupInterpolatorProvider("Lookup");

//////////////////////////////////////////////////////////////////////////////

void LookupInterpolator::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< string >
    ("InputFileName","Input data file from which all variables will be extrapolated");
}
    
//////////////////////////////////////////////////////////////////////////////

LookupInterpolator::LookupInterpolator(const std::string& name) :
  StateInterpolator(name)
{
  addConfigOptionsTo(this);

  m_infile = "";
  setParameter("InputFileName",&m_infile); 
}


//////////////////////////////////////////////////////////////////////////////

LookupInterpolator::~LookupInterpolator()
{
  cf_assert(m_lookupState.size() == 0);
  cf_assert(m_lookupState.capacity() == 0);
}

//////////////////////////////////////////////////////////////////////////////

void LookupInterpolator::setup()
{
  CFAUTOTRACE;
  
  const std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  
  StateInterpolator::setup(); 
  
  if (m_infile != "") {
    if (PE::GetPE().IsParallel()) {
      PE::GetPE().setBarrier(nsp);
      for (CFuint p = 0; p < PE::GetPE().GetProcessorCount(nsp); ++p) {
	if (p == PE::GetPE().GetRank (nsp)) {
	  fillTable();
	}
	PE::GetPE().setBarrier(nsp);
      }
    }
    else {
      fillTable();
    }
  }
  else {
    CFLog(WARN, "WARNING: LookupInterpolator::setup() => filename not specified!\n");
  } 
}

//////////////////////////////////////////////////////////////////////////////

void LookupInterpolator::unsetup()
{
  CFAUTOTRACE;

  StateInterpolator::unsetup(); 
  
  for (CFuint i = 0; i < m_lookupState.size(); ++i) {
    deletePtr(m_lookupState[i]);
  }
  SwapEmpty(m_lookupState);
}
    
//////////////////////////////////////////////////////////////////////////////

void LookupInterpolator::fillTable()
{
  boost::filesystem::path filepath = Environment::DirPaths::getInstance().
    getWorkingDir() / m_infile;
  Common::SelfRegistPtr<Environment::FileHandlerInput>* fhandle =
    Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().createPtr();
  ifstream& fin = (*fhandle)->open(filepath);
  
  // this interpolator expects a file with the following format
  //
  // N (number of point to read)
  // s1 v11 v12 ... v1Neq
  // s2 v21 v22 ... v2Neq
  // ...
  // sN vN1 vN2 ... vNNeq
  //
  // where "s" is a (curvilinear) coordinate and "Neq" is the number of equations
  
  string variables;
  // read the first line with the variables names
  getline(fin,variables);
  CFuint nbPoints = 0;
  fin >> nbPoints;
  
  // allocate the look up tables
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  m_lookupState.resize(nbEqs);
  for (CFuint i = 0; i < nbEqs; ++i) {
    m_lookupState[i] = new LookUpTable<CFreal, CFreal>(nbPoints);
  } 
  
  // nbEqs + "y" AL: hardcoded here
  CFreal ycoord = 0.;
  CFreal tmpVar = 0.;
  for (CFuint ip = 0; ip < nbPoints; ++ip) {
    fin >> ycoord;
    CFLog(DEBUG_MAX, "LookupInterpolator::fillTable() => " << ycoord << " ");
    for (CFuint i = 0; i < nbEqs; ++i) {
      fin >> tmpVar;
      m_lookupState[i]->insert(ycoord, tmpVar);
      CFLog(DEBUG_MAX, tmpVar << " ");
    }
    CFLog(DEBUG_MAX, "\n");
  }
  
  // sort the data 
  for (CFuint i = 0; i < nbEqs; ++i) {
    m_lookupState[i]->sortKeys();
  }
  
  fin.close();
  delete fhandle;
}

//////////////////////////////////////////////////////////////////////////////

void LookupInterpolator::interpolate(const CFuint varID, 
				     const CFreal& in, CFreal& out)
{
  out = m_lookupState[varID]->get(in);
}
    
//////////////////////////////////////////////////////////////////////////////

void LookupInterpolator::interpolate(const RealVector& in, RealVector& out)
{
  throw NotImplementedException(FromHere(),"LookupInterpolator::interpolate() missing");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
