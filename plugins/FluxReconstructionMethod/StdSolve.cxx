// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"

#include "FluxReconstructionMethod/StdSolve.hh"
#include "FluxReconstructionMethod/FluxReconstruction.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

namespace COOLFluiD {
  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StdSolve,FluxReconstructionSolverData,FluxReconstructionModule >
  stdSolveProvider("StdSolve");
  
//////////////////////////////////////////////////////////////////////////////
  
StdSolve::StdSolve(const std::string& name) :
  FluxReconstructionSolverCom(name),
  socket_gradients("gradients"),
  socket_normals("normals")
  {
  }
  
//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSink > >
StdSolve::needsSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSink > > result;
  result.push_back(&socket_gradients);
  result.push_back(&socket_normals);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdSolve::execute()
{
  CFAUTOTRACE;
  
  CFLog(INFO, "StdSolve::execute()\n");
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod
}  // namespace COOLFluiD

