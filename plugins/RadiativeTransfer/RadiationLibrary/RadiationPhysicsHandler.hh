#ifndef COOLFluiD_RadiativeTransfer_RadiationPhysicsHandler_hh
#define COOLFluiD_RadiativeTransfer_RadiationPhysicsHandler_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/OwnedObject.hh"
#include "Common/SetupObject.hh"
#include "Common/NonCopyable.hh"
#include "Environment/ConcreteProvider.hh"
#include "RadiativeTransfer/RadiationLibrary/RadiationPhysics.hh"
#include "Framework/MethodCommand.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

namespace RadiativeTransfer {

class RadiationPhysics;

class RadiationPhysicsHandler : public Common::OwnedObject,
                                public Config::ConfigObject,
                                public Common::NonCopyable<RadiationPhysicsHandler>
{
public:

  typedef Environment::ConcreteProvider<RadiationPhysicsHandler,1> PROVIDER;
  typedef const std::string& ARG1;


  RadiationPhysicsHandler(const std::string &name);
  ~RadiationPhysicsHandler();

  static std::string getClassName() { return "RadiationPhysicsHandler"; }
  void configure(Config::ConfigArgs& args);
  static void defineConfigOptions(Config::OptionList& options);

  void setupWavStride(CFuint loop);

  Common::SharedPtr< RadiationPhysics > getCellDistPtr(CFuint stateID);
  Common::SharedPtr< RadiationPhysics > getWallDistPtr(CFuint GhostStateID);


  CFuint getNumberLoops() const {return m_nbLoops;}

  void setupDataSockets(Framework::SocketBundle sockets){
      m_sockets = sockets;
  }

  Framework::SocketBundle* const getDataSockets(){return &m_sockets;}

  void setup();
    
  void configureTRS();
  CFuint getNbStates();
  void setupAxiFlag(bool isAxi){ m_isAxi = isAxi;}

  inline void getWallTRSnames(std::vector<std::string> &wallTRSnames)
         {wallTRSnames=m_wallTRSnames;}

  inline void getBoundaryTRSnames(std::vector<std::string> &boundaryTRSnames)
         {boundaryTRSnames = m_boundaryTRSnames;}

  inline void getMediumTRSnames(std::vector<std::string> &mediumTRSnames)
          {mediumTRSnames = m_mediumTRSnames;}

  inline bool isStateNull(CFuint stateID) {
    return (m_statesOwner[stateID][0] == -1);
  }

  inline bool isGhostStateNull(CFuint ghostStateID) {
    return (m_ghostStatesOwner[ghostStateID][0] == -1);
  }


  inline CFuint getCurrentCellStateID(){return m_cellStateID;}
  inline CFuint getCurrentCellTrsIdx(){return m_cellStateOwnerIdx;}

  inline CFuint getCurrentWallGhostStateID(){return m_ghostStateID;}
  inline CFuint getCurrentWallTrsIdx(){return m_ghostStateOwnerIdx;}
  inline CFuint getCurrentWallGeoID(){return m_ghostStateWallGeoID;}

  inline CFuint getTempID(){return m_TempID;}
  inline CFuint getNbTemps(){return m_nbTemps;}
  inline bool isAxi(){return m_isAxi;}

  CFuint getNbGhostStates();
private:
  Framework::SocketBundle m_sockets;
  std::vector< std::string > m_radiationPhysicsNames;
  std::vector<Common::SharedPtr< RadiationPhysics > > m_radiationPhysics;
  std::vector<std::vector<CFint> > m_statesOwner, m_ghostStatesOwner;

  std::vector<std::string> m_wallTRSnames;
  std::vector<std::string> m_boundaryTRSnames;
  std::vector<std::string> m_mediumTRSnames;

  CFreal m_wavMin;
  CFreal m_wavMax;
  CFuint m_nbLoops;
  CFint m_TempID;
  CFuint m_nbTemps;

  CFuint m_cellStateID;
  CFuint m_cellStateOwnerIdx;

  CFuint m_ghostStateID;
  CFuint m_ghostStateOwnerIdx;
  CFuint m_ghostStateWallGeoID;

  bool m_isAxi;

};

}
}

#endif
