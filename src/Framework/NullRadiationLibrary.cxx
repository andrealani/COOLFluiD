// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/NullRadiationLibrary.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/Framework.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  
  namespace Framework {
        
//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NullRadiationLibrary,
			    RadiationLibrary,
			    FrameworkLib,
			    1>
nullRadiationLibraryProvider("Null");

//////////////////////////////////////////////////////////////////////////////

NullRadiationLibrary::NullRadiationLibrary(const std::string& name) :
  RadiationLibrary(name)
{
}
    
//////////////////////////////////////////////////////////////////////////////
    
NullRadiationLibrary::~NullRadiationLibrary()
{
}
 
//////////////////////////////////////////////////////////////////////////////
    
void NullRadiationLibrary::setup()
{
}
    
//////////////////////////////////////////////////////////////////////////////
    
void NullRadiationLibrary::runOnStagnationLine(Common::SafePtr<std::vector<CFuint> > stagnationLineCells,
					       Framework::ProxyDofIterator<CFreal>* pstates,
					       CFreal* qrad)
{
}

//////////////////////////////////////////////////////////////////////////////
    
void NullRadiationLibrary::runOnStructuredMesh(const std::vector<std::vector<CFuint>* >& meshByLine,
					       Framework::ProxyDofIterator<CFreal>* pstates,
					       CFreal* qrad)
{
}

//////////////////////////////////////////////////////////////////////////////

void NullRadiationLibrary::computeProperties(Framework::ProxyDofIterator<CFreal>* pstates,
					     RealMatrix& data, CFuint iWavRange)
{
}

/////////////////////////////////////////////////////////////////////////////
 
} // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
