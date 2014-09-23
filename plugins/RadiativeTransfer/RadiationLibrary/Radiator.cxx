#include "Radiator.hh"

namespace COOLFluiD {

namespace RadiativeTransfer {

/// Constructor without arguments
Radiator::Radiator(const std::string& name):
        Common::OwnedObject(),
        ConfigObject(name)
{
  addConfigOptionsTo(this);
}


CFreal Radiator::getCellVolume( CFuint stateID )
{
      static Framework::DataHandle<CFreal> volumes
          = m_radPhysicsHandlerPtr->getDataSockets()->volumes.getDataHandle();
      static Framework::DataHandle<Framework::State*, Framework::GLOBAL> states
              = m_radPhysicsHandlerPtr->getDataSockets()->states.getDataHandle();
      static bool isAxi = m_radPhysicsHandlerPtr->isAxi();
      //static CFuint trsTypeID  = m_radPhysicsPtr->getTRStypeID();
      CFreal volume = volumes[stateID];
      //convert to axi volumes if necessary
      volume *= (isAxi)? 2.*3.141516*(*states[stateID]).getCoordinates()[YY] : 1. ;
      return volume;
}


CFreal Radiator::getCurrentCellVolume()
{
      CFuint stateID = m_radPhysicsHandlerPtr->getCurrentCellStateID();
      return getCellVolume( stateID );
}


CFreal Radiator::getWallArea( CFuint wallGeoID)
{
      //Framework::State* state = m_radPhysicsPtr->getWallState( wallTrsIdx );

      static Framework::DataHandle<CFreal> faceAreas
          = m_radPhysicsHandlerPtr->getDataSockets()->faceAreas.getDataHandle();

      static Framework::DataHandle<CFreal> faceCenters
          = m_radPhysicsHandlerPtr->getDataSockets()->faceCenters.getDataHandle();

      static bool isAxi = m_radPhysicsHandlerPtr->isAxi();
      //convert to axi volumes if necessary
      //std::cout<<"area: "<<faceAreas[wallGeoID]<<std::endl;
      return faceAreas[wallGeoID] * ( (isAxi)? 2.*3.141516* faceCenters[wallGeoID*DIM_2D+YY] : 1. );
}

CFreal Radiator::getCurrentWallArea()
{
  CFuint wallGeoID = m_radPhysicsHandlerPtr->getCurrentWallGeoID();
  //CFuint wallTrsIdx = m_radPhysicsHandlerPtr->getCurrentWallTrsIdx();

  return getWallArea(wallGeoID);
}

}

}
