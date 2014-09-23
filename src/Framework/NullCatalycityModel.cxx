// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "NullCatalycityModel.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {
    
//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<NullCatalycityModel, CatalycityModel, FrameworkLib, 1>
nullCatalycityModelProv("Null");

//////////////////////////////////////////////////////////////////////////////

NullCatalycityModel::NullCatalycityModel(const std::string& name)
  : CatalycityModel(name)
{
}
    
//////////////////////////////////////////////////////////////////////////////

NullCatalycityModel::~NullCatalycityModel()
{
}

//////////////////////////////////////////////////////////////////////////////

void NullCatalycityModel::computeMassProduction(const CFreal temp, 
						const CFreal rho,
						const RealVector& ys, 
						RealVector& mp)
{
  CFout << "NullCatalycityModel::computeMassProduction() called!" << "\n";  
}

//////////////////////////////////////////////////////////////////////////////
  
} // namespace Framework
  
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
