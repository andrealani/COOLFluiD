#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"
#include "Common/BadValueException.hh"

#include "StructMech/Materials.hh"
#include "StructMech/MaterialData.hh"
#include "StructMech/StructMechTerm.hh"
#include "StructMech/StructMechPhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

void StructMechPhysicalModel::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("MaterialNames","Material Names.");
  options.addConfigOption< std::vector<std::string> >("ZoneNames","Name of the various zones of the mesh.");
}

//////////////////////////////////////////////////////////////////////////////

StructMechPhysicalModel::StructMechPhysicalModel(const std::string& name)
  : PhysicalModelImpl(name),
    m_zonesMap(),
    m_zoneID(0),
    m_structTerm(new StructMechTerm("StructTerm"))
{
  addConfigOptionsTo(this);

  m_materialStr = std::vector<std::string>();
  setParameter("MaterialNames",&m_materialStr);

  m_zoneNames = std::vector<std::string>();
  setParameter("ZoneNames",&m_zoneNames);


}

//////////////////////////////////////////////////////////////////////////////

StructMechPhysicalModel::~StructMechPhysicalModel()
{
}

//////////////////////////////////////////////////////////////////////////////

void StructMechPhysicalModel::configure ( Config::ConfigArgs& args )
{
  PhysicalModelImpl::configure(args);

  if(m_zoneNames.size() == 0){
    m_zoneNames.resize(1);
    m_zoneNames[0] = "InnerCells";
  }
  if(m_materialStr.size() == 0){
    m_materialStr.resize(1);
    m_materialStr[0] = "Custom";
  }

  cf_assert(m_materialStr.size() == m_zoneNames.size());
  m_material.resize(m_materialStr.size());

  for(CFuint iMat=0; iMat < m_material.size(); iMat++)
  {
    m_material[iMat] =
      Environment::Factory<MaterialPropertyLib>::getInstance().
      getProvider(m_materialStr[iMat])->create(m_materialStr[iMat]);

    configureNested(m_material[iMat].getPtr(), args);
    cf_assert(m_material[iMat].isNotNull());
  }

}

//////////////////////////////////////////////////////////////////////////////

void StructMechPhysicalModel::setup()
{
  PhysicalModelImpl::setup();

  for(CFuint iZone=0; iZone < m_zoneNames.size(); ++iZone){
    m_zonesMap.insert(m_zoneNames[iZone], iZone);
  }
  m_zonesMap.sortKeys();

}

//////////////////////////////////////////////////////////////////////////////

RealMatrix* StructMechPhysicalModel::getStiffnessMatrix()
{
  return m_material[m_zoneID]->computeStiffnessMatrix();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Framework::BaseTerm>
StructMechPhysicalModel::getConvectiveTerm() const
{
  return m_structTerm.get();
}

//////////////////////////////////////////////////////////////////////////////

Common::SafePtr<Framework::BaseTerm>
StructMechPhysicalModel::getDiffusiveTerm() const
{
  return m_structTerm.get();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
