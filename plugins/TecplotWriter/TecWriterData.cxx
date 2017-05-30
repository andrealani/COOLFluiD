// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "TecplotWriter/TecplotWriter.hh"
#include "TecplotWriter/TecWriterData.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/CFLog.hh"
#include "Framework/ConvectiveVarSet.hh"
#include "Framework/VarSetTransformer.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Framework/SpaceMethod.hh"
#include "Framework/SpaceMethodData.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace TecplotWriter {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<TecWriterData>, TecWriterData, TecplotWriterModule> nullTecWriterComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void TecWriterData::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< bool >("printExtraValues","flag telling if to print extra values");
   options.addConfigOption< std::string >("outputVar","output variables set");
   options.addConfigOption< std::vector<std::string> >("SurfaceTRS","List of TRS's to be writen in the surface file.");
   options.addConfigOption< bool >("SurfaceOnly","Print only the surface data chosen in SurfaceTRS"); 
   options.addConfigOption< bool >("CoordinatesOnly","Print only the coordinates");
   options.addConfigOption< bool >("WithEquations","Print the equations");
   options.addConfigOption< bool >("AppendAuxData","Boolean switch to append 'AUXDATA' fields in zone headers or not.");
}
      
//////////////////////////////////////////////////////////////////////////////

TecWriterData::TecWriterData(Common::SafePtr<Framework::Method> owner)
  : OutputFormatterData(owner),
    m_filepath(),
    m_outputVarStr(),
    m_outputVarSet()
{
  addConfigOptionsTo(this);
  
  m_outputVarStr = "Null";
  setParameter("outputVar",&m_outputVarStr);

  m_printExtraValues = false;
  setParameter("printExtraValues",&m_printExtraValues);

  m_surface_only = false;
  setParameter("SurfaceOnly",&m_surface_only);
  
  m_coord_only = false;
  setParameter("CoordinatesOnly",&m_coord_only);
  
  m_withEquations = true;
  setParameter("WithEquations",&m_withEquations);
  
  m_surfTRS = std::vector<std::string>();
  setParameter("SurfaceTRS",&m_surfTRS);

  m_appendAuxData = false;
  setParameter("AppendAuxData",&m_appendAuxData);
}

//////////////////////////////////////////////////////////////////////////////

TecWriterData::~TecWriterData()
{
}

//////////////////////////////////////////////////////////////////////////////

void TecWriterData::configure ( Config::ConfigArgs& args )
{
  OutputFormatterData::configure(args);
  
  // configure the update varset -------------------------------------------
  
  std::string updateVarStr = "Null";
  SafePtr<SpaceMethodData> smData = getCollaborator<SpaceMethod>()->getSpaceMethodData();
  if (smData.isNotNull()) {
    updateVarStr   = smData->getUpdateVarStr();
    m_updateVarSet = smData->getUpdateVar();
  }
  
  if (m_outputVarStr == "Null") { m_outputVarStr = updateVarStr; }
  
  // get the physical model -------------------------------------------
  
  std::string name = getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);
  
  std::string provider_base = (physModel->getConvectiveName() != "Null") ? physModel->getConvectiveName() : physModel->getDiffusiveName() ;
  
  // configure the output varset -------------------------------------------

  CFLog(VERBOSE, "TecWriterData::outputVarStr = " << provider_base + m_outputVarStr << "\n");
  
  const std::string varProv = provider_base + m_outputVarStr;
  m_outputVarSet.reset
    (FACTORY_GET_PROVIDER(getFactoryRegistry(), ConvectiveVarSet, varProv)->
     create(physModel->getImplementor()->getConvectiveTerm()));
  
  cf_assert(m_outputVarSet.isNotNull());
  
  if (updateVarStr == "Null") {
    updateVarStr = m_outputVarStr; 
    m_updateVarSet = m_outputVarSet.getPtr();
  }
  
  cf_assert ( updateVarStr != "Null" && m_outputVarStr != "Null" );
  
  // create the transformer from update to output variables ----------------
  std::string provider_tr =  Framework::VarSetTransformer::getProviderName
    (physModel->getConvectiveName(), updateVarStr, m_outputVarStr);
  m_updateToOutputVar = 
    FACTORY_GET_PROVIDER(getFactoryRegistry(), VarSetTransformer, provider_tr)->
    create(physModel->getImplementor());
  
  cf_assert(m_updateToOutputVar.isNotNull());  
}
      
//////////////////////////////////////////////////////////////////////////////

void TecWriterData::setup ()
{
  CFAUTOTRACE;
  
  OutputFormatterData::setup();
  
  m_updateToOutputVar->setup(1);
  m_outputVarSet->setup();
}
      
//////////////////////////////////////////////////////////////////////////////

void TecWriterData::unsetup ()
{
  CFAUTOTRACE;
  
  OutputFormatterData::unsetup();
}

//////////////////////////////////////////////////////////////////////////////
      
    } // namespace TecplotWriter

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

