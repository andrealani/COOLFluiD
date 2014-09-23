// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/MethodCommandProvider.hh"
#include "CFmeshFileWriter/CFmeshWriterData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace CFmeshFileWriter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<CFmeshWriterData>, CFmeshWriterData, CFmeshFileWriterModule> nullCFmeshWriterComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void CFmeshWriterData::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< bool >("StorePastStates","Flag to store the 'pastStates' values.");
   options.addConfigOption< bool >("StorePastNodes","Flag to store the 'pastNodes' values.");
   options.addConfigOption< bool >("StoreInterStates","Flag to store the 'interStates' values.");
   options.addConfigOption< bool >("StoreInterNodes","Flag to store the 'interNodes' values.");
   options.addConfigOption< std::vector<std::string> >("ExtraNodalVarNames","Extra Nodal Data to write to the file.");
   options.addConfigOption< std::vector<std::string> >("ExtraStateVarNames","Extra State Data to write to the file.");
   options.addConfigOption< std::vector<std::string> >("ExtraVarNames","Extra Data to write to the file.");
   options.addConfigOption< std::vector<CFuint> >("ExtraNodalVarStrides","Extra Nodal Data Strides to write to the file.");
   options.addConfigOption< std::vector<CFuint> >("ExtraStateVarStrides","Extra State Data Strides to write to the file.");
   options.addConfigOption< std::vector<CFuint> >("ExtraVarStrides","Extra Data Strides to write to the file.");

   //   options.addConfigOption< std::vector<std::string> >("ExtraNodalVarTags","Tag in CFmesh for extra Nodal Data to read from the file.");
   //   options.addConfigOption< std::vector<std::string> >("ExtraStateVarTags","Tag in CFmesh for extra State Data to read from the file.");
}

//////////////////////////////////////////////////////////////////////////////

CFmeshWriterData::CFmeshWriterData(Common::SafePtr<Framework::Method> owner)
 : OutputFormatterData(owner)
{
   addConfigOptionsTo(this);

   _storePastStates = false;
   setParameter("StorePastStates", &_storePastStates);

   _storePastNodes = false;
   setParameter("StorePastNodes", &_storePastNodes);
   
  _storeInterStates = false;
   setParameter("StoreInterStates", &_storeInterStates);

   _storeInterNodes = false;
   setParameter("StoreInterNodes", &_storeInterNodes);

   _extraNodalVarNames = std::vector<std::string>();
   setParameter("ExtraNodalVarNames", &_extraNodalVarNames);

   _extraStateVarNames = std::vector<std::string>();
   setParameter("ExtraStateVarNames", &_extraStateVarNames);

   _extraNodalVarStrides = std::vector<CFuint>();
   setParameter("ExtraNodalVarStrides", &_extraNodalVarStrides);

   _extraStateVarStrides = std::vector<CFuint>();
   setParameter("ExtraStateVarStrides", &_extraStateVarStrides);
   
   _extraVarNames = std::vector<std::string>();
   setParameter("ExtraVarNames",&_extraVarNames);
   
   _extraVarStrides = std::vector<CFuint>();
   setParameter("ExtraVarStrides", &_extraVarStrides);

   ///@todo for the moment the tags are not used
   //    _extraNodalVarTags = std::vector<std::string>();
   //    setParameter("ExtraNodalVarTags", &_extraNodalVarTags);
   //
   //    _extraStateVarTags = std::vector<std::string>();
   //    setParameter("ExtraStateVarTags", &_extraStateVarTags);
}

//////////////////////////////////////////////////////////////////////////////

CFmeshWriterData::~CFmeshWriterData()
{
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshWriterData::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  OutputFormatterData::configure(args);

  ///Check the size of the vectors for the extra variables to read from the file
  if(_extraNodalVarNames.size() > 0)
  {
    if(_extraNodalVarTags.size() > 0){
      cf_assert(_extraNodalVarTags.size() == _extraNodalVarNames.size());
      cf_assert(_extraNodalVarStrides.size() == _extraNodalVarNames.size());
    }
    else{
      cf_assert(_extraNodalVarStrides.size() == _extraNodalVarNames.size());
      _extraNodalVarTags.resize(_extraNodalVarNames.size());
      for(CFuint iVar = 0; iVar < _extraNodalVarNames.size(); iVar++){
        _extraNodalVarTags[iVar] = _extraNodalVarNames[iVar];
      }
    }
  }

  if(_extraStateVarNames.size() > 0)
  {
    if(_extraStateVarTags.size() > 0){
      cf_assert(_extraStateVarTags.size() == _extraStateVarNames.size());
      cf_assert(_extraStateVarStrides.size() == _extraStateVarNames.size());
    }
    else{
      cf_assert(_extraStateVarStrides.size() == _extraStateVarNames.size());
      _extraStateVarTags.resize(_extraStateVarNames.size());
      for(CFuint iVar = 0; iVar < _extraStateVarNames.size(); iVar++){
        _extraStateVarTags[iVar] = _extraStateVarNames[iVar];
      }
    }
  }
  
  if(_extraVarNames.size() > 0)
  {
    if(_extraVarTags.size() > 0){
      cf_assert(_extraVarTags.size() == _extraVarNames.size());
      cf_assert(_extraVarStrides.size() == _extraVarNames.size());
    }
    else{
      cf_assert(_extraVarStrides.size() == _extraVarNames.size());
      _extraVarTags.resize(_extraVarNames.size());
      for(CFuint iVar = 0; iVar < _extraVarNames.size(); iVar++){
        _extraVarTags[iVar] = _extraVarNames[iVar];
      }
    }
  }

  // Build map between the File TAG and the socket NAME
  for(CFuint iVar = 0; iVar < _extraNodalVarNames.size(); iVar++){
    _extraNodalVarMap.insert(_extraNodalVarTags[iVar],pair<std::string,CFuint>
			     (_extraNodalVarNames[iVar],_extraNodalVarStrides[iVar]));
  }

  for(CFuint iVar = 0; iVar < _extraStateVarNames.size(); iVar++){
    _extraStateVarMap.insert(_extraStateVarTags[iVar],pair<std::string,CFuint>
			     (_extraStateVarNames[iVar],_extraStateVarStrides[iVar]));
  }
  
  for(CFuint iVar = 0; iVar < _extraVarNames.size(); iVar++){
    _extraVarMap.insert(_extraVarTags[iVar],pair<std::string,CFuint>
			     (_extraVarNames[iVar],_extraVarStrides[iVar]));
  }

  _extraNodalVarMap.sortKeys();
  _extraStateVarMap.sortKeys();
  _extraVarMap.sortKeys();
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshWriterData::setup()
{
  OutputFormatterData::setup();
}

//////////////////////////////////////////////////////////////////////////////

void CFmeshWriterData::unsetup()
{
  OutputFormatterData::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace CFmeshFileWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

