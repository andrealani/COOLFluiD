// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>
#include <iomanip>

#include "Common/COOLFluiD.hh"
#include "Common/PE.hh"

#include "Common/CFMap.hh"
#include "Environment/FileHandlerOutput.hh"

#include "Common/OSystem.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/DirPaths.hh"

#include "Framework/MapGeoEnt.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/DataHandleOutput.hh"
#include "Framework/SubSystemStatus.hh"

#include "TecplotWriter/TecplotWriter.hh"
#include "TecplotWriter/WriteSolution1D.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Common;
using namespace boost::filesystem;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WriteSolution1D, TecWriterData, TecplotWriterModule>
writeSolution1DProvider("WriteSolution1D");

//////////////////////////////////////////////////////////////////////////////

WriteSolution1D::WriteSolution1D(const std::string& name) :
  WriteSolution(name)
{
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolution1D::execute()
{
  CFLog(INFO, "Writing solution to: " << getMethodData().getFilename().string() << "\n");
  writeToFile(getMethodData().getFilename());
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolution1D::writeToFileStream(std::ofstream& fout)
{
  CFAUTOTRACE;

  SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<ProxyDofIterator<RealVector>*> nstatesProxy =
    socket_nstatesProxy.getDataHandle();

  // this isa sort of handle for the nodal states
  // (which can be stored as arrays of State*, RealVector* or
  // RealVector but they are used as arrays of RealVector*)
  ProxyDofIterator<RealVector>& nodalStates = *nstatesProxy[0];

  // we will assume that the number of nodes is the same as
  // the number of states but the connectivity might be different
  //   cf_assert(nodes.size() == nodalStates.getSize());

  SafePtr<DataHandleOutput> datahandle_output = getMethodData().getDataHOutput();
  datahandle_output->getDataHandles();

// unused //  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal refL = PhysicalModelStack::getActive()->getImplementor()->getRefLength();

  //  Tecplot Header
  fout << "TITLE      =  \"Unstructured grid data\"" << "\n";
  fout << "VARIABLES  = \"x\" ";

  SafePtr<ConvectiveVarSet> updateVarSet = getMethodData().getUpdateVarSet();
  const vector<std::string>& varNames = updateVarSet->getVarNames();
  cf_assert(varNames.size() == nbEqs);

  for (CFuint i = 0 ;  i < nbEqs; ++i) {
    std::string n = varNames[i];
    if ( *n.begin()   != '\"' )  n  = '\"' + n;
    if ( *(n.end()--) != '\"' )  n += '\"';
    fout << " " << n;
  }

  if (getMethodData().shouldPrintExtraValues()) {
    vector<std::string> extraVarNames = updateVarSet->getExtraVarNames();
    for (CFuint i = 0 ;  i < extraVarNames.size(); ++i) {
      fout << " " << extraVarNames[i];
    }
  }

  datahandle_output->printVarNames(fout);

  fout << "\n";

  // array to store the dimensional states
  RealVector dimState(nbEqs);
  RealVector extraValues; // size will be set in the VarSet
  State tempState;

  const CFuint nbLocalNodes = nodes.size();
  for (CFuint iNode = 0; iNode < nbLocalNodes; ++iNode) {
    const CFuint nodeID = nodes[iNode]->getLocalID();

    // node has to be printed with the right length
    fout << setw(30) << fixed << setprecision(30) << (*nodes[iNode])[XX]*refL << " ";

    const RealVector& currState = *nodalStates.getState(nodeID);
    for (CFuint ieq = 0; ieq < nbEqs; ++ieq) {
      tempState[ieq] = currState[ieq];
    }

    const CFuint stateID = nodalStates.getStateLocalID(nodeID);
    tempState.setLocalID(stateID);
    // the node is set  in the temporary state
    tempState.setSpaceCoordinates(nodes[nodeID]);

    if (getMethodData().shouldPrintExtraValues()) {
      // dimensionalize the solution
      updateVarSet->setDimensionalValuesPlusExtraValues
	(tempState, dimState, extraValues);
      fout << dimState << " " << extraValues << "\n";
    }
    else {
      // set other useful (dimensional) physical quantities
      updateVarSet->setDimensionalValues(tempState, dimState);
      fout << dimState << "\n";
    }
  }
  fout.close();
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolution1D::setup()
{
  CFAUTOTRACE;

  WriteSolution::setup();
}

//////////////////////////////////////////////////////////////////////////////


    } // namespace TecplotWriter

} // namespace COOLFluiD
