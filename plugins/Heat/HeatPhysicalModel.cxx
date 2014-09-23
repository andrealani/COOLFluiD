#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"
#include "Common/BadValueException.hh"

#include "Heat/HeatPhysicalModel.hh"
#include "Heat/ConvHeatTerm.hh"
#include "Heat/DiffHeatTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Heat {

//////////////////////////////////////////////////////////////////////////////

void HeatPhysicalModel::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Conductivity","Heat conductivity parameter.");
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
  options.addConfigOption< std::vector<std::string> >("ZoneNames","Name of the various zones of the mesh.");
}

//////////////////////////////////////////////////////////////////////////////

HeatPhysicalModel::HeatPhysicalModel(const std::string& name)
  : PhysicalModelImpl(name),
    m_physicalData(),
    m_convectiveTerm(new ConvHeatTerm("ConvTerm")),
    m_diffusiveTerm(new DiffHeatTerm("DiffTerm")),
    m_constantConductivity(true),
    m_Computedconductivity(),
    m_zonesMap(),
    m_zoneID(0)
{
   addConfigOptionsTo(this);
   m_conductivity = 1.0;
   setParameter("Conductivity",&m_conductivity);

   m_functions = std::vector<std::string>();
   setParameter("Def",&m_functions);

   m_vars = std::vector<std::string>();
   setParameter("Vars",&m_vars);

   m_zoneNames = std::vector<std::string>();
   setParameter("ZoneNames",&m_zoneNames);

}

//////////////////////////////////////////////////////////////////////////////

void HeatPhysicalModel::configure ( Config::ConfigArgs& args )
{
  PhysicalModelImpl::configure(args);

  if((m_functions.size() > 0) && (m_vars.size() > 0)){
    m_constantConductivity = false;
    m_vFunction.setFunctions(m_functions);
    m_vFunction.setVariables(m_vars);
    try {
      m_vFunction.parse();
    }
    catch (Common::ParserException& e) {
      CFout << e.what() << "\n";
      throw; // retrow the exception to signal the error to the user
    }
  }

  if ( m_conductivity <= 0.0)
    throw Common::BadValueException (FromHere(),"Conductivity must be greater than zero.");
}

//////////////////////////////////////////////////////////////////////////////

HeatPhysicalModel::~HeatPhysicalModel()
{
}

//////////////////////////////////////////////////////////////////////////////

void HeatPhysicalModel::setup()
{
  PhysicalModelImpl::setup();

  if(m_zoneNames.size() == 0){
    m_zoneNames.resize(1);
    m_zoneNames[0] = "InnerCells";
  }
  if(!m_constantConductivity) cf_assert(m_zoneNames.size() == m_functions.size());

  for(CFuint iZone=0; iZone < m_zoneNames.size(); ++iZone){
    m_zonesMap.insert(m_zoneNames[iZone], iZone);
  }
  m_zonesMap.sortKeys();
  m_Computedconductivity.resize(m_zoneNames.size());
}

//////////////////////////////////////////////////////////////////////////////

void HeatPhysicalModel::setCurrentZone(const std::string zoneName)
{
  PhysicalModelImpl::setCurrentZone(zoneName);

  m_zoneID = m_zonesMap.find(zoneName);

}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Framework::BaseTerm>
HeatPhysicalModel::getConvectiveTerm() const
{
  return m_convectiveTerm.get();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Framework::BaseTerm>
HeatPhysicalModel::getDiffusiveTerm() const
{
  return m_diffusiveTerm.get();
}

//////////////////////////////////////////////////////////////////////////////

CFreal HeatPhysicalModel::getConductivity(const RealVector& coord, const RealVector& state)
{

  if(!m_constantConductivity){
    //Evaluate the function
    for(CFuint iDim=0; iDim < m_nbDim; ++iDim) m_variables[iDim] = coord[iDim];
    for(CFuint iEq=0; iEq < m_nbEqs; ++iEq) m_variables[m_nbDim + iEq] = state[iEq];

    m_vFunction.evaluate(m_variables,m_Computedconductivity);
/*CFout << "m_Computedconductivity: " << m_Computedconductivity << "\n";
CFout << "m_Computedconductivity["<<m_zoneID<<": " << m_Computedconductivity[m_zoneID] << "\n";*/
    return m_Computedconductivity[m_zoneID];
  }

  return m_conductivity;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Heat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
