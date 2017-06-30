// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MeshFormatConverter.hh"
#include "Common/Stopwatch.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/SingleBehaviorFactory.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

void MeshFormatConverter::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("Discontinuous","Is the solution space be discontinuous?");
  options.addConfigOption< vector<std::string> >("IgnoreTRS",
	"Names of the TRS to ignore while converting");
}

//////////////////////////////////////////////////////////////////////////////

MeshFormatConverter::MeshFormatConverter(const std::string& name)
 :  Common::OwnedObject(),
    ConfigObject(name),
    Common::NullableObject()
{
  addConfigOptionsTo(this);

  m_isDiscontinuous = false;
  setParameter("Discontinuous",&m_isDiscontinuous);

  m_ignoreTRSNames = vector<std::string>();
  setParameter("IgnoreTRS",&m_ignoreTRSNames);
}

//////////////////////////////////////////////////////////////////////////////

MeshFormatConverter::~MeshFormatConverter()
{
}

//////////////////////////////////////////////////////////////////////////////

void MeshFormatConverter::configure ( Config::ConfigArgs& args )
{
  ConfigObject::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void MeshFormatConverter::convert(const boost::filesystem::path& fromFilepath,
				  const boost::filesystem::path& filepath)
{
  CFAUTOTRACE;

  Stopwatch<WallTime> stp;
  stp.start();

  // only reads if not yet been read
  readFiles(fromFilepath);
  
  stp.stop();
  CFLog(INFO, "Reading " << this->getName() << " took: " << stp.read() << "s\n");
  
  stp.start();
  adjustToCFmeshNodeNumbering();

  Common::SelfRegistPtr<Environment::FileHandlerOutput>* fhandle =
    Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().createPtr();
  ofstream& fout = (*fhandle)->open(filepath);
  
  // AL: Those info make the CFmesh file version-dependent, can create problems with regression testing 
  fout << "!COOLFLUID_VERSION "    << Environment::CFEnv::getInstance().getCFVersion() << "\n";
  // this can fail if there are problems with SVN
  // fout << "!COOLFLUID_SVNVERSION " << Environment::CFEnv::getInstance().getSvnVersion() << "\n";
  fout << "!CFMESH_FORMAT_VERSION 1.3\n";
  
  fout << "!NB_DIM " << getDimension() << "\n";

  CFuint nbVariables = getNbVariables();
  if (nbVariables == 0) {
    nbVariables = PhysicalModelStack::getActive()->getNbEq();
  }

  fout << "!NB_EQ " << nbVariables << "\n";

  if (isDiscontinuous())
  {
    writeDiscontinuousElements(fout);
    writeDiscontinuousTrsData(fout);
    writeNodes(fout);
    writeDiscontinuousStates(fout);
  }
  else
  {
    writeContinuousElements(fout);
    writeContinuousTrsData(fout);
    writeNodes(fout);
    writeContinuousStates(fout);
  }
  fout << "!END" << "\n";

  stp.stop();
  (*fhandle)->close();

  delete fhandle;
  CFLog(INFO, "Conversion " << this->getName()<< " took: " << stp.read() << "s\n");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
