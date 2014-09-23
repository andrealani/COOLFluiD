#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"
#include "Common/BadValueException.hh"

#include "StructMechHeat/StructMechHeatPhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::StructMech;
//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMechHeat {

//////////////////////////////////////////////////////////////////////////////

void StructMechHeatPhysicalModel::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Conductivity","Heat conductivity parameter.");
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
  options.addConfigOption< CFreal >("InitialTemp","Initial Temp such that u=v=0.");
}

//////////////////////////////////////////////////////////////////////////////

StructMechHeatPhysicalModel::StructMechHeatPhysicalModel(const std::string& name)
  : StructMechPhysicalModel(name),
    m_constantConductivity(true),
    m_Computedconductivity(1)
{
   addConfigOptionsTo(this);
   m_conductivity = 1.0;
   setParameter("Conductivity",&m_conductivity);

   m_functions = std::vector<std::string>();
   setParameter("Def",&m_functions);

   m_vars = std::vector<std::string>();
   setParameter("Vars",&m_vars);

   m_initialTemp = 298.15;
   setParameter("InitialTemp",&m_initialTemp);

}

//////////////////////////////////////////////////////////////////////////////

StructMechHeatPhysicalModel::~StructMechHeatPhysicalModel()
{
}

//////////////////////////////////////////////////////////////////////////////

void StructMechHeatPhysicalModel::configure ( Config::ConfigArgs& args )
{
  StructMechPhysicalModel::configure(args);

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

  if(m_conductivity <= 0.0) {
    throw Common::BadValueException (FromHere(),"Conductivity must be greater than zero.");
  }

}

//////////////////////////////////////////////////////////////////////////////

void StructMechHeatPhysicalModel::setup()
{
  StructMechPhysicalModel::setup();

  m_variables.resize(m_nbDim + m_nbEqs);

  if(!m_constantConductivity) cf_assert(m_zoneNames.size() == m_functions.size());
  m_Computedconductivity.resize(m_zoneNames.size());

}

//////////////////////////////////////////////////////////////////////////////

void StructMechHeatPhysicalModel::setCurrentZone(const std::string zoneName)
{
  StructMechPhysicalModel::setCurrentZone(zoneName);

  m_zoneID = m_zonesMap.find(zoneName);

}

//////////////////////////////////////////////////////////////////////////////

CFreal StructMechHeatPhysicalModel::getConductivity(const RealVector& coord, const RealVector& state)
{
  if(!m_constantConductivity){
    //Evaluate the function
    for(CFuint iDim=0; iDim < m_nbDim; ++iDim) m_variables[iDim] = coord[iDim];
    for(CFuint iEq=0; iEq < m_nbEqs; ++iEq) m_variables[m_nbDim + iEq] = state[iEq];

    m_vFunction.evaluate(m_variables,m_Computedconductivity);

    return m_Computedconductivity[m_zoneID];
  }

  return m_conductivity;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMechHeat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
