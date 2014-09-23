// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "MarcoTest/MarcoTest.hh"
#include "MarcoTest/MarcoTestMethodData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace MarcoTest {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<MarcoTestMethodData>, 
		      MarcoTestMethodData, 
		      MarcoTestModule> 
nullMarcoTestMethodComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void MarcoTestMethodData::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >("GeoDataComputer","Computer of geometric data (e.g. normals, volumes).");
  options.addConfigOption< std::string >("NodalExtrapolation","Nodal extrapolation strategy.");
}

//////////////////////////////////////////////////////////////////////////////

MarcoTestMethodData::MarcoTestMethodData(Common::SafePtr<Framework::Method> owner) : 
  ConvergenceMethodData(owner),
  m_faceTrsGeoBuilder(),
  m_faceCellTrsGeoBuilder(),
  m_cellTrsGeoBuilder(),
  m_geoWithNodesBuilder(),
  m_geoDataComputer(),
  m_nStatesExtrapolator()
{
  addConfigOptionsTo(this);
  
  m_geoDataComputerStr = "FVMCC";
  setParameter("GeoDataComputer",&m_geoDataComputerStr);
  
  m_nStatesExtrapolatorStr = "DistanceBased";
  setParameter("NodalExtrapolation", &m_nStatesExtrapolatorStr);
}
      
//////////////////////////////////////////////////////////////////////////////

MarcoTestMethodData::~MarcoTestMethodData()
{
  CFLog(VERBOSE, "MarcoTestMethodData::~MarcoTestMethodData()\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void MarcoTestMethodData::configure ( Config::ConfigArgs& args )
{
  ConvergenceMethodData::configure(args);

  // configure the geometric data computer
  configureGeoDataComputer(args);
  
  // configure the nodal extrapolator 
  configureNodalStatesExtrapolator(args);
}

//////////////////////////////////////////////////////////////////////////////

void MarcoTestMethodData::setup()
{
  ConvergenceMethodData::setup();
  
  m_faceTrsGeoBuilder.setup();
  m_faceCellTrsGeoBuilder.setup();
  m_cellTrsGeoBuilder.setup();
  m_geoWithNodesBuilder.setup();
}
      
//////////////////////////////////////////////////////////////////////////////

void MarcoTestMethodData::unsetup()
{
  m_faceTrsGeoBuilder.unsetup();
  m_faceCellTrsGeoBuilder.unsetup();
  m_cellTrsGeoBuilder.unsetup();
  m_geoWithNodesBuilder.unsetup();
  
  ConvergenceMethodData::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void MarcoTestMethodData::configureGeoDataComputer ( Config::ConfigArgs& args )
{
  std::string name = m_geoDataComputerStr;
  
  SharedPtr<MarcoTestMethodData> thisPtr(this);
  Common::SafePtr<BaseMethodStrategyProvider<MarcoTestMethodData,GeoDataComputer<MarcoTestMethodData> > > prov =
    Environment::Factory<GeoDataComputer<MarcoTestMethodData> >::getInstance().getProvider(name);
  
  cf_assert(prov.isNotNull());
  
  m_geoDataComputer = prov->create(name,thisPtr);
  configureNested ( m_geoDataComputer.getPtr(), args ); 
}

//////////////////////////////////////////////////////////////////////////////

void MarcoTestMethodData::configureNodalStatesExtrapolator ( Config::ConfigArgs& args )
{
  // create the nodal states extrapolator
  SharedPtr<MarcoTestMethodData> thisPtr(this);
  
  Common::SafePtr<BaseMethodStrategyProvider<
    MarcoTestMethodData, NodalStatesExtrapolator<MarcoTestMethodData> > > prov =
    Environment::Factory<NodalStatesExtrapolator<MarcoTestMethodData> >::getInstance().
    getProvider(m_nStatesExtrapolatorStr);
  
  cf_assert(prov.isNotNull());
  m_nStatesExtrapolator = prov->create(m_nStatesExtrapolatorStr,thisPtr);
  configureNested ( m_nStatesExtrapolator.getPtr(), args );
  cf_assert(m_nStatesExtrapolator.isNotNull());
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MarcoTest

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

