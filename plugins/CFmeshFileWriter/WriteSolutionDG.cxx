// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "CFmeshFileWriter/CFmeshFileWriter.hh"
#include "Common/SafePtr.hh"
#include "WriteSolutionDG.hh"
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

MethodCommandProvider<WriteSolutionDG, CFmeshWriterData, CFmeshFileWriterModule> WriteSolutionDGprovider("WriteSolutionDG");

//////////////////////////////////////////////////////////////////////////////

WriteSolutionDG::WriteSolutionDG(const std::string& name) :
  CFmeshWriterCom(name),
  m_writer(),
  m_data(new CFmeshWriterSource()),
  m_sockets()
{
  SafePtr<CFmeshWriterSource> ptr = m_data.get();
  m_writer.setWriteData(ptr);
}

//////////////////////////////////////////////////////////////////////////////

WriteSolutionDG::~WriteSolutionDG()
{
}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionDG::configure ( Config::ConfigArgs& args )
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
}

//////////////////////////////////////////////////////////////////////

void WriteSolutionDG::setup()
{
  CFAUTOTRACE;

  m_data->setMeshData();
  m_data->consistencyCheck();
  m_data->setPastDataStorageFlags(getMethodData().storePastNodes(),
                                  getMethodData().storePastStates() );
  m_data->setExtraVarNamesAndTags(getMethodData().getExtraNVarSocketNamesAndTags(),
                                  getMethodData().getExtraSVarSocketNamesAndTags(),
                                  getMethodData().getExtraVarSocketNamesAndTags());

  m_data->setExtraDataSockets(&m_sockets);


}

//////////////////////////////////////////////////////////////////////////////

void WriteSolutionDG::execute()
{
  CFAUTOTRACE;

  // Write the File
  m_writer.writeToFile(getMethodData().getFilename());
  
  CFLog(INFO, "Writing solution to: " << getMethodData().getFilename().string() << "\n");
}
      
//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
WriteSolutionDG::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = m_sockets.getAllSinkSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileWriter

} // namespace COOLFluiD
