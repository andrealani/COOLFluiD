// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "CFmeshFileWriter/CFmeshFileWriter.hh"
#include "CFmeshFileWriter/ParWriteSolution.hh"
#include "CFmeshFileWriter/ParCFmeshFileWriter.hh"
#include "CFmeshFileWriter/ParCFmeshBinaryFileWriter.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileWriter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ParWriteSolution<ParCFmeshFileWriter>, 
		      CFmeshWriterData, CFmeshFileWriterModule>
parWriteSolutionProvider("ParWriteSolution");

MethodCommandProvider<ParWriteSolution<ParCFmeshBinaryFileWriter>, 
		      CFmeshWriterData, CFmeshFileWriterModule>
parWriteBinarySolutionProvider("ParWriteBinarySolution");

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileWriter

} // namespace COOLFluiD
