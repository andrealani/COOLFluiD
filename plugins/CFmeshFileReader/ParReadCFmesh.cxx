// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Environment/ConcreteProvider.hh"
#include "Framework/MethodCommandProvider.hh"
#include "CFmeshFileReader/CFmeshFileReader.hh"
#include "CFmeshFileReader/ParReadCFmesh.hh"
#include "CFmeshFileReader/ParCFmeshFileReader.hh"
#include "CFmeshFileReader/ParCFmeshBinaryFileReader.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileReader {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ParReadCFmesh<ParCFmeshFileReader>, 
		      CFmeshReaderData, 
		      CFmeshFileReaderPlugin>
stdParReadCFmeshProvider("ParReadCFmesh");

MethodCommandProvider<ParReadCFmesh<ParCFmeshBinaryFileReader>, 
		      CFmeshReaderData, 
		      CFmeshFileReaderPlugin>
stdParReadCFmeshBinaryProvider("ParReadCFmeshBinary");

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileReader

} // namespace COOLFluiD
