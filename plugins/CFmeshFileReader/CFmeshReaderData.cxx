// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "CFmeshFileReader/CFmeshFileReader.hh"
#include "CFmeshFileReader/CFmeshReaderData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileReader {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<CFmeshReaderData>,
		      CFmeshReaderData,
		      CFmeshFileReaderPlugin>
nullCFmeshReaderComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderData::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("convertFromFile","Name of file in other format to convert From.");
   options.addConfigOption< std::vector<CFreal> >("TranslationVector","Vector to translate the mesh.");
   options.addConfigOption< bool >("TranslateMesh","Flag to say if the mesh is to be translated.");

   options.addConfigOption< bool >("AnisotropicScaling","Flag to say if the mesh is to be scaled anisotropicaly.");
   options.addConfigOption< std::vector<CFreal> >("ScalingVector","Vector to scale the mesh.");
   options.addConfigOption< std::vector<bool> >
     ("UseInitValues",
      "Specifies which variables will be initialized with the given values.");

   options.addConfigOption< std::vector<CFreal> >
    ("InitValues",
     "Initial values for the specified variables.");

   options.addConfigOption< std::vector<CFuint> >
    ("InitValuesIDs",
     "Initial values IDs for the specified variables.");

   options.addConfigOption< bool >("ReadPastStates","Flag to read the 'pastStates' values.");
   options.addConfigOption< bool >("ReadPastNodes","Flag to read the 'pastNodes' values.");
   options.addConfigOption< bool >("ReadInterStates","Flag to read the 'interStates' values.");
   options.addConfigOption< bool >("ReadInterNodes","Flag to read the 'interNodes' values.");
   options.addConfigOption< std::vector<std::string> >("ExtraNodalVarNames","Extra Nodal Data to read from the file.");
   options.addConfigOption< std::vector<std::string> >("ExtraStateVarNames","Extra State Data to read from the file.");
   options.addConfigOption< std::vector<std::string> >("ExtraVarNames","Extra Data to read from the file.");
   options.addConfigOption< std::vector<std::string> >("ExtraNodalVarTags","Tag in CFmesh for extra Nodal Data to read from the file.");
   options.addConfigOption< std::vector<std::string> >("ExtraStateVarTags","Tag in CFmesh for extra State Data to read from the file.");
   options.addConfigOption< std::vector<std::string> >("ExtraVarTags","Tag in CFmesh for extra Data to read from the file.");
   options.addConfigOption< std::vector<CFuint> >("ExtraNodalVarStrides","Extra Nodal Data Strides to read from the file.");
   options.addConfigOption< std::vector<CFuint> >("ExtraStateVarStrides","Extra State Data Strides to read from the file.");
   options.addConfigOption< std::vector<CFuint> >("ExtraVarStrides","Extra Data Strides to read from the file.");
   options.addConfigOption< std::vector<CFreal> >("StateScalingValues","Values to scale the state variables.");
   options.addConfigOption< std::string >("SolutionFile","Name of solution file from which each processor can read itw own states."); 
   options.addConfigOption< std::vector<CFuint> >("NodeSwitchIDs","IDs for switching the node coordinates");
}

//////////////////////////////////////////////////////////////////////////////

CFmeshReaderData::CFmeshReaderData(Common::SafePtr<Framework::Method> owner)
  : MeshCreatorData(owner)
{
  addConfigOptionsTo(this);

  m_convertFromFileName = "";
  setParameter("convertFromFile",&m_convertFromFileName);

  m_isTranslated = false;
  setParameter("TranslateMesh",&m_isTranslated);

  m_translationVector = std::vector<CFreal>();
  setParameter("TranslationVector",&m_translationVector);
   m_isAnisotropicScaling = false;
   setParameter("AnisotropicScaling",&m_isAnisotropicScaling);

   m_scalingVector = std::vector<CFreal>();
   setParameter("ScalingVector",&m_scalingVector);

  m_useInitValues = std::vector<bool>();
  setParameter("UseInitValues", &m_useInitValues);

  m_initValues = std::vector<CFreal>();
  setParameter("InitValues", &m_initValues);

  m_initValuesIDs = std::vector<CFuint>();
  setParameter("InitValuesIDs", &m_initValuesIDs);

  m_readPastStates = false;
  setParameter("ReadPastStates", &m_readPastStates);

  m_readPastNodes = false;
  setParameter("ReadPastNodes", &m_readPastNodes);

  m_readInterStates = false;
  setParameter("ReadInterStates", &m_readInterStates);

  m_readInterNodes = false;
  setParameter("ReadInterNodes", &m_readInterNodes);

  m_extraNodalVarNames = std::vector<std::string>();
  setParameter("ExtraNodalVarNames", &m_extraNodalVarNames);

  m_extraStateVarNames = std::vector<std::string>();
  setParameter("ExtraStateVarNames", &m_extraStateVarNames);

  m_extraVarNames = std::vector<std::string>();
  setParameter("ExtraVarNames", &m_extraVarNames);

  m_extraNodalVarTags = std::vector<std::string>();
  setParameter("ExtraNodalVarTags", &m_extraNodalVarTags);

  m_extraStateVarTags = std::vector<std::string>();
  setParameter("ExtraStateVarTags", &m_extraStateVarTags);

  m_extraVarTags = std::vector<std::string>();
  setParameter("ExtraVarTags", &m_extraVarTags);

  m_extraNodalVarStrides = std::vector<CFuint>();
  setParameter("ExtraNodalVarStrides", &m_extraNodalVarStrides);

  m_extraStateVarStrides = std::vector<CFuint>();
  setParameter("ExtraStateVarStrides", &m_extraStateVarStrides);

  m_extraVarStrides = std::vector<CFuint>();
  setParameter("ExtraVarStrides", &m_extraVarStrides);

  m_stateScalingValues = std::vector<CFreal>();
  setParameter("StateScalingValues", &m_stateScalingValues);

   m_solutionFile = "Null";
  setParameter("SolutionFile",&m_solutionFile);

  m_nodeSwitchIDs = std::vector<CFuint>();
  setParameter("NodeSwitchIDs", &m_nodeSwitchIDs);
}

//////////////////////////////////////////////////////////////////////////////

CFmeshReaderData::~CFmeshReaderData()
{
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderData::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  MeshCreatorData::configure(args);

  // if not other file name was set to convert from,
  // tale as defaut the one supplied for reading.
  if(m_convertFromFileName == "")
  {
    m_convertFromFileName = m_fileName;
  }
  else
  {
    cf_assert(m_extraNodalVarNames.size() == 0);
    cf_assert(m_extraStateVarNames.size() == 0);
    cf_assert(m_extraVarNames.size() == 0);
  }

  // Check the size of the vectors for the extra variables to read from the file
  if(m_extraNodalVarNames.size() > 0)
  {
    if(m_extraNodalVarTags.size() > 0){
      cf_assert(m_extraNodalVarTags.size() == m_extraNodalVarNames.size());
      cf_assert(m_extraNodalVarStrides.size() == m_extraNodalVarNames.size());
    }
    else{
      cf_assert(m_extraNodalVarStrides.size() == m_extraNodalVarNames.size());
      m_extraNodalVarTags.resize(m_extraNodalVarNames.size());
      for(CFuint iVar = 0; iVar < m_extraNodalVarNames.size(); iVar++){
        m_extraNodalVarTags[iVar] = m_extraNodalVarNames[iVar];
      }
    }
  }

  if(m_extraStateVarNames.size() > 0)
  {
    if(m_extraStateVarTags.size() > 0){
      cf_assert(m_extraStateVarTags.size() == m_extraStateVarNames.size());
      cf_assert(m_extraStateVarStrides.size() == m_extraStateVarNames.size());
    }
    else{
      cf_assert(m_extraStateVarStrides.size() == m_extraStateVarNames.size());
      m_extraStateVarTags.resize(m_extraStateVarNames.size());
      for(CFuint iVar = 0; iVar < m_extraStateVarNames.size(); iVar++){
        m_extraStateVarTags[iVar] = m_extraStateVarNames[iVar];
      }
    }
  }
  
  if(m_extraVarNames.size() > 0)
  {
    if(m_extraVarTags.size() > 0){
      cf_assert(m_extraVarTags.size() == m_extraVarNames.size());
      cf_assert(m_extraVarStrides.size() == m_extraVarNames.size());
    }
    else{
      cf_assert(m_extraVarStrides.size() == m_extraVarNames.size());
      m_extraVarTags.resize(m_extraVarNames.size());
      for(CFuint iVar = 0; iVar < m_extraVarNames.size(); iVar++){
        m_extraVarTags[iVar] = m_extraVarNames[iVar];
      }
    }
  }

  // Build map between the File TAG and the socket NAME
  for(CFuint iVar = 0; iVar < m_extraNodalVarNames.size(); iVar++){
    cf_assert(m_extraNodalVarTags.size() == m_extraNodalVarNames.size());
    cf_assert(m_extraNodalVarStrides.size() == m_extraNodalVarNames.size());

    m_extraNodalVarMap.insert(m_extraNodalVarTags[iVar],pair<std::string,CFuint>
			      (m_extraNodalVarNames[iVar],m_extraNodalVarStrides[iVar]));
  }

  for(CFuint iVar = 0; iVar < m_extraStateVarNames.size(); iVar++){
    cf_assert(m_extraStateVarTags.size() == m_extraStateVarNames.size());
    cf_assert(m_extraStateVarStrides.size() == m_extraStateVarNames.size());
    m_extraStateVarMap.insert(m_extraStateVarTags[iVar],pair<std::string,CFuint>
			     (m_extraStateVarNames[iVar],m_extraStateVarStrides[iVar]));
  }
  
  for(CFuint iVar = 0; iVar < m_extraVarNames.size(); iVar++){
    cf_assert(m_extraVarTags.size() == m_extraVarNames.size());
    cf_assert(m_extraVarStrides.size() == m_extraVarNames.size());
    m_extraVarMap.insert(m_extraVarTags[iVar],pair<std::string,CFuint>
			     (m_extraVarNames[iVar],m_extraVarStrides[iVar]));
  }
  
  m_extraNodalVarMap.sortKeys();
  m_extraStateVarMap.sortKeys();
  m_extraVarMap.sortKeys();

}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderData::setup()
{
  MeshCreatorData::setup();
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshReaderData::unsetup()
{
  MeshCreatorData::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileReader

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

