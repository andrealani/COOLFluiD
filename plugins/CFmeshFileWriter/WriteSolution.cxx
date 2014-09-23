// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "CFmeshFileWriter/CFmeshFileWriter.hh"
#include "Common/SafePtr.hh"
#include "WriteSolution.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/DirPaths.hh"
#include "Common/FilesystemException.hh"
#include "Framework/MethodCommandProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileWriter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WriteSolution, CFmeshWriterData, CFmeshFileWriterModule> writeSolutionProvider("WriteSolution");

//////////////////////////////////////////////////////////////////////////////

WriteSolution::WriteSolution(const std::string& name) :
  CFmeshWriterCom(name),
  m_writer(),
  m_data(new CFmeshWriterSource()),
  m_sockets()
{
  SafePtr<CFmeshWriterSource> ptr = m_data.get();
  m_writer.setWriteData(ptr);
}

//////////////////////////////////////////////////////////////////////////////

WriteSolution::~WriteSolution()
{
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolution::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  CFmeshWriterCom::configure(args);

  Common::SafePtr<Common::CFMap<std::string, pair<std::string,CFuint> > > extraNVars =
    getMethodData().getExtraNVarSocketNamesAndTags();
  Common::SafePtr<Common::CFMap<std::string, pair<std::string,CFuint> > > extraSVars =
    getMethodData().getExtraSVarSocketNamesAndTags();

  for(CFuint iNVar = 0; iNVar < extraNVars->size();iNVar++){
    std::string socketName = (*extraNVars)[iNVar].first;
    m_sockets.createSocketSink<CFreal>(socketName);
  }

  for(CFuint iSVar = 0; iSVar < extraSVars->size();iSVar++){
    std::string socketName = (*extraSVars)[iSVar].first;
    m_sockets.createSocketSink<CFreal>(socketName);
  }

  bool storePastNodes = getMethodData().storePastNodes();
  bool storePastStates = getMethodData().storePastStates();

  if(storePastNodes) m_sockets.createSocketSink<Node*>("pastNodes");
  if(storePastStates) m_sockets.createSocketSink<State*>("pastStates");

  bool storeInterNodes = getMethodData().storeInterNodes();
  bool storeInterStates = getMethodData().storeInterStates();

  if(storeInterNodes) m_sockets.createSocketSink<Node*>("interNodes");
  if(storeInterStates) m_sockets.createSocketSink<State*>("interStates");

}
//////////////////////////////////////////////////////////////////////

void WriteSolution::setup()
{
  CFAUTOTRACE;

  m_data->setMeshData();
  m_data->consistencyCheck();
  m_data->setPastDataStorageFlags(getMethodData().storePastNodes(),
                                  getMethodData().storePastStates() );
  
  m_data->setInterDataStorageFlags(getMethodData().storeInterNodes(),
                                  getMethodData().storeInterStates() );

  m_data->setExtraVarNamesAndTags(getMethodData().getExtraNVarSocketNamesAndTags(),
                                  getMethodData().getExtraSVarSocketNamesAndTags(),
                                  getMethodData().getExtraVarSocketNamesAndTags());

  m_data->setExtraDataSockets(&m_sockets);

}

//////////////////////////////////////////////////////////////////////////////

void WriteSolution::execute()
{
  CFAUTOTRACE;

  // Write the File
  m_writer.writeToFile(getMethodData().getFilename());

  CFout << "Writing solution to: " << getMethodData().getFilename().string() << "\n";
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WriteSolution::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = m_sockets.getAllSinkSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileWriter

} // namespace COOLFluiD
