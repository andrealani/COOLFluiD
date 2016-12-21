#include "RadiativeTransfer/RadiationLibrary/Radiator.hh"
#include "RadiativeTransfer/RadiationLibrary/RadiationPhysicsHandler.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

void Radiator::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint >("nbBands", "Number of Bands");
  options.addConfigOption< CFuint >("nbBins", "Number of Bins");
}

//////////////////////////////////////////////////////////////////////////////

/// Constructor without arguments
Radiator::Radiator(const std::string& name):
  Common::OwnedObject(),
  ConfigObject(name),
  m_angstrom(1e-10),
  m_states(CFNULL),
  m_volumes(CFNULL),
  m_faceAreas(CFNULL),
  m_faceCenters(CFNULL)
{
  addConfigOptionsTo(this);
  
  m_nbBands = 1;
  setParameter("nbBands", &m_nbBands);
  
  m_nbBins = 1;
  setParameter("nbBins", &m_nbBins);
}

//////////////////////////////////////////////////////////////////////////////

void Radiator::setup()
{
  m_states = m_radPhysicsHandlerPtr->getDataSockets()->states.getDataHandle();
  m_volumes = m_radPhysicsHandlerPtr->getDataSockets()->volumes.getDataHandle();

  m_faceAreas = m_radPhysicsHandlerPtr->getDataSockets()->faceAreas.getDataHandle();
  
  m_faceCenters = m_radPhysicsHandlerPtr->getDataSockets()->faceCenters.getDataHandle();
}

//////////////////////////////////////////////////////////////////////////////

CFreal Radiator::getCurrentCellVolume() const
{
  const CFuint stateID = 
    m_radPhysicsHandlerPtr->getCurrentCellStateID();
  return getCellVolume( stateID );
}

//////////////////////////////////////////////////////////////////////////////

CFreal Radiator::getCurrentWallArea() const 
{
  const CFuint wallGeoID = 
      m_radPhysicsHandlerPtr->getCurrentWallGeoID();
  //CFuint wallTrsIdx = m_radPhysicsHandlerPtr->getCurrentWallTrsIdx();
  return getWallArea(wallGeoID);
}
  
//////////////////////////////////////////////////////////////////////////////
  
CFreal Radiator::getCellVolume(CFuint stateID) const
{
  CFreal volume = m_volumes[stateID];
  cf_assert(volume > 0.);
  if (m_radPhysicsHandlerPtr->isAxi()) { //convert to axi volumes if necessary
    volume *= 2.*MathTools::MathConsts::CFrealPi()*(*m_states[stateID]).getCoordinates()[YY];
  }
  return volume;
}

//////////////////////////////////////////////////////////////////////////////
  
CFreal Radiator::getWallArea(CFuint wallGeoID) const
{  
  CFreal area = m_faceAreas[wallGeoID];
  cf_assert(area > 0.);
  if (m_radPhysicsHandlerPtr->isAxi()) { //convert to axi volumes if necessary
    const CFreal h = m_faceCenters[wallGeoID*DIM_2D+YY];
    const CFreal coeff = h * 2.*MathTools::MathConsts::CFrealPi();
    area *= coeff;
  }
  //std::cout<<" r: "<<h<<" coeff: "<<coeff<<" final area: "<<area<<std::endl;
  return area ;
}

//////////////////////////////////////////////////////////////////////////////

}

}
