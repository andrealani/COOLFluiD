// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "SimpleGlobalMeshAdapter/SimpleGlobalMeshAdapter.hh"
#include "CopyFilesPrepare.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Environment/DirPaths.hh"
#include "Common/OSystem.hh"
#include "Common/FilesystemException.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace boost::filesystem;

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SimpleGlobalMeshAdapter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<CopyFilesPrepare, SimpleMeshAdapterData, SimpleGlobalMeshAdapterModule> copyFilesPrepareProvider("CopyFilesPrepare");

//////////////////////////////////////////////////////////////////////////////

void CopyFilesPrepare::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("InitialFiles","Name of the files to copy from working to result dir.");
}

//////////////////////////////////////////////////////////////////////////////

CopyFilesPrepare::CopyFilesPrepare(const std::string& name)  :
  SimpleMeshAdapterCom(name)
{

   addConfigOptionsTo(this);

   _initialFilesStr = std::vector<std::string>();
   setParameter("InitialFiles",&_initialFilesStr);

}

//////////////////////////////////////////////////////////////////////////////

void CopyFilesPrepare::execute()
{

  //Here you want to copy files from working directory to the result dir
  ///@todo do this only at initialization
  if(true)
  {
    cf_assert(_initialFilesStr.size() > 0);

    // Copying the initial files to the working dir
    for(CFuint iFile=0;iFile < _initialFilesStr.size();iFile++)
    {
      path user_file = path(_initialFilesStr[iFile]);

      boost::filesystem::path from_file = DirPaths::getInstance().getWorkingDir() / user_file;
      boost::filesystem::path to_file   = DirPaths::getInstance().getResultsDir() / user_file.filename();

      Common::SelfRegistPtr<Environment::FileHandlerInput> fileHandle =
        Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();

      // this forces the download of the file if it does not exist
      fileHandle->open(from_file);
      fileHandle->close();

      // delete destination file if exsits to avoid exception throw
      if ( boost::filesystem::exists (to_file) )
        boost::filesystem::remove (to_file);

      boost::filesystem::copy_file( from_file, to_file );
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace SimpleGlobalMeshAdapter

  } // namespace Numerics

} // namespace COOLFluiD
