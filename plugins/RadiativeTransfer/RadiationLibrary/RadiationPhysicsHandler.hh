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

  /// get the number of loops
  CFuint getNumberLoops() const {return m_nbLoops;}

  /// set up the data sockets
  void setupDataSockets(Framework::SocketBundle sockets){m_sockets = sockets;}
  
  /// get the bundle of sockets
  Framework::SocketBundle* const getDataSockets(){return &m_sockets;}
  
  /// setup private data
  void setup();
  
  /// configure the corresponing TRS
  void configureTRS();
  
  /// get the number of states
  CFuint getNbStates() const {return m_statesOwner.size();}
  
  /// set up the axisymmetric flag
  void setupAxiFlag(const bool isAxi)
  {
    if (isAxi) {
      cf_assert(Framework::PhysicalModelStack::getActive()->getDim() == DIM_2D);
    }
    m_isAxi = isAxi;
  }
  
  /// get the Wall TRS names
  void getWallTRSnames(std::vector<std::string> &wallTRSnames)
  {wallTRSnames=m_wallTRSnames;}
  
  /// get the Boundary TRS names
  void getBoundaryTRSnames(std::vector<std::string> &boundaryTRSnames)
  {boundaryTRSnames = m_boundaryTRSnames;}
  
  /// get the Medium TRS names
  void getMediumTRSnames(std::vector<std::string> &mediumTRSnames)
  {mediumTRSnames = m_mediumTRSnames;}
  
  /// check if the given state is NULL
  bool isStateNull(CFuint stateID) const 
  {
    return (m_statesOwner[stateID][0] == -1);
  }
  
  /// check if the given ghost state is NULL
  bool isGhostStateNull(CFuint ghostStateID) const 
  {
    return (m_ghostStatesOwner[ghostStateID][0] == -1);
  }
  
  /// @return the current cell state ID
  CFuint getCurrentCellStateID() const {return m_cellStateID;}
  
  /// @return the current cell ID into its corresponding  TRS
  CFuint getCurrentCellTrsIdx() const {return m_cellStateOwnerIdx;}
  
  /// @return the current cell wall ghost state ID
  CFuint getCurrentWallGhostStateID() const {return m_ghostStateID;}
  
  /// @return the current wall face ID into its corresponding TRS
  CFuint getCurrentWallTrsIdx() const {return m_ghostStateOwnerIdx;}
  
  /// @return the current wall face ID (local ID in the current processor)
  CFuint getCurrentWallGeoID() const {return m_ghostStateWallGeoID;}
  
  /// @return the variable ID corresponding to the temperature
  CFuint getTempID() const {return m_TempID;}
  
  /// @return the numer of temperature
  /// @post   will return an integer>1 for thermal nonequilibrium models
  CFuint getNbTemps() const {return m_nbTemps;}

  /// flag telling whether the case is 2D axialsymemtric
  bool isAxi() const {return m_isAxi;}
  
  /// @return the number of ghost states
  CFuint getNbGhostStates() const {return m_ghostStatesOwner.size();} 
  
private:
  Framework::SocketBundle m_sockets;
  std::vector<Common::SharedPtr< RadiationPhysics > > m_radiationPhysics;
  std::vector<std::vector<CFint> > m_statesOwner;
  std::vector<std::vector<CFint> > m_ghostStatesOwner;

  std::vector<std::string> m_wallTRSnames;
  std::vector<std::string> m_boundaryTRSnames;
  std::vector<std::string> m_mediumTRSnames;
  
  CFuint m_nbTemps;
  CFuint m_cellStateID;
  CFuint m_cellStateOwnerIdx;
  
  CFuint m_ghostStateID;
  CFuint m_ghostStateOwnerIdx;
  CFuint m_ghostStateWallGeoID;
  
  bool m_isAxi;
  
  CFreal m_wavMin;
  CFreal m_wavMax;
  CFuint m_nbLoops;
  CFint  m_TempID;
  
  std::vector< std::string > m_radiationPhysicsNames;
};

//////////////////////////////////////////////////////////////////////////////

}
}

#endif
