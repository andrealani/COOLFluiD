// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"

#include "CFmeshFileReader/CFmeshFileReader.hh"
#include "CFmeshFileReader/ReadDummy.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileReader {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ReadDummy, CFmeshReaderData, CFmeshFileReaderPlugin>
readDummyCFmeshProvider("Dummy");

//////////////////////////////////////////////////////////////////////////////

ReadDummy::ReadDummy(const std::string& name) :
  CFmeshReaderCom(name),
  m_sockets(),
  socket_nodes("nodes"),
  socket_states("states")
{
}

//////////////////////////////////////////////////////////////////////////////

ReadDummy::~ReadDummy()
{
}

//////////////////////////////////////////////////////////////////////////////

void ReadDummy::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);
}

/////////////////////////////////////////////////////////////////////////////

void ReadDummy::execute()
{
  DataHandle<State*,GLOBAL> states = socket_states.getDataHandle(); 
  states.reserve(1, PhysicalModelStack::getActive()->getNbEq()*sizeof(CFreal));
  DataHandle<Node*,GLOBAL> nodes = socket_nodes.getDataHandle();
  nodes.reserve(1, PhysicalModelStack::getActive()->getDim()*sizeof(CFreal));
}
      
//////////////////////////////////////////////////////////////////////////////
      
std::vector<Common::SafePtr<BaseDataSocketSink> >
ReadDummy::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = m_sockets.getAllSinkSockets();

  result.push_back(&socket_nodes);
  result.push_back(&socket_states);

  return result;
}

/////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileReader

} // namespace COOLFluiD
