// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/FileHandlerInputConcrete.hh"
#include "Environment/ObjectProvider.hh"
#include "Environment/DirectFileAccess.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace Environment {

//////////////////////////////////////////////////////////////////////////////

// provider for this behavior
Environment::ObjectProvider<
                Environment::FileHandlerInputConcrete< DirectFileAccess >,
                Environment::FileHandlerInput,
                EnvironmentModule >
directFileAccessProvider("DirectFileAccess");

//////////////////////////////////////////////////////////////////////////////

    } // namespace Environment

} // namespace COOLFluiD
