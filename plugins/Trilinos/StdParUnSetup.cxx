// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Trilinos/Trilinos.hh"
#include "StdParUnSetup.hh"
#include "Framework/MeshData.hh"
#include "MathTools/RealVector.hh"
#include "Framework/MethodCommandProvider.hh"

#include "Epetra_Comm.h"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace Trilinos {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdParUnSetup, TrilinosLSSData, TrilinosModule> stdParUnSetupProvider("StdParUnSetup");

//////////////////////////////////////////////////////////////////////////////

void StdParUnSetup::execute()
{
  CFAUTOTRACE;

  // map
  Epetra_Map *map = getMethodData().getEpetraMap();
  if (map!=NULL) {delete map; map=NULL;}

  // communicator
  Epetra_Comm *comm = getMethodData().getEpetraComm();
  if (comm!=NULL) {delete comm; comm=NULL;}
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Trilinos

  } // namespace Numerics

} // namespace COOLFluiD
