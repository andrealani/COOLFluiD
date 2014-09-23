#include "Environment/DirPaths.hh"

#include "Framework/BadFormatException.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "NavierStokes/Euler3DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

#include "SpectralFV/BaseBndFaceTermComputer.hh"
#include "SpectralFV/SpectralFVElementData.hh"

#include "SpectralFVNavierStokes/BCCurvedWallEuler3D.hh"
#include "SpectralFVNavierStokes/SpectralFVNavierStokes.hh"

#include "Common/ConnectivityTable.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider< BCCurvedWallEuler3D,
                                   SpectralFVMethodData,
                                   BCStateComputer,
                                   SpectralFVNavierStokesModule
                                 > BCCurvedWallEuler3DProvider("CurvedWallEuler3D");

//////////////////////////////////////////////////////////////////////////////

void BCCurvedWallEuler3D::defineConfigOptions(Config::OptionList& options)
{
  CF_DEBUG_POINT;
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("NormalsDef","Definition of the Functions that describe the normals.");
}

//////////////////////////////////////////////////////////////////////////////

BCCurvedWallEuler3D::BCCurvedWallEuler3D(const std::string& name) :
  BCStateComputer(name),
  m_nodes(CFNULL),
  m_eulerVarSet(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData(),
  m_flxPntWheightCoordsSVFaces(CFNULL),
  m_faceFlxPntNormals()
{
  addConfigOptionsTo(this);

  m_functions = vector<std::string>();
  setParameter("NormalsDef",&m_functions);

  m_vars = vector<std::string>();
  setParameter("Vars",&m_vars);
}

//////////////////////////////////////////////////////////////////////////////

BCCurvedWallEuler3D::~BCCurvedWallEuler3D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCCurvedWallEuler3D::computeGhostStates(const vector< State* >& intStates,
                                             vector< State* >& ghostStates,
                                             const std::vector< RealVector >& normals,
                                             const std::vector< RealVector >& coords)
{
  // number of states
  const CFuint nbrStates = ghostStates.size();
  cf_assert(nbrStates == intStates.size());
  cf_assert(nbrStates == normals.size());

  // get some physical data from the model
  const CFreal gamma = m_eulerVarSet->getModel()->getGamma();
  const CFreal gammaDivGammaMinus1 = gamma/(gamma-1.0);

  // get the real normals for the current face
  vector< RealVector >& flxPntNormals = m_faceFlxPntNormals[m_faceID];

  // loop over the states
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    // normal
    const RealVector& normal = flxPntNormals[iState];

    // dereference states
    State& intState   = (*intStates[iState]);
    State& ghostState = (*ghostStates[iState]);

    cf_assert(intState.size() == 5);
    cf_assert(ghostState.size() == 5);

    // set the physical data starting from the inner state
    m_eulerVarSet->computePhysicalData(intState,m_intSolPhysData);

    // compute normal velocity component
    const CFreal uNX2 = 2.0*(m_intSolPhysData[EulerTerm::VX]*normal[XX] +
                             m_intSolPhysData[EulerTerm::VY]*normal[YY] +
                             m_intSolPhysData[EulerTerm::VZ]*normal[ZZ]);

    // set the physical data for the ghost state
    m_ghostSolPhysData[EulerTerm::RHO] = m_intSolPhysData[EulerTerm::RHO];
    m_ghostSolPhysData[EulerTerm::VX]  = m_intSolPhysData[EulerTerm::VX] - uNX2*normal[XX];
    m_ghostSolPhysData[EulerTerm::VY]  = m_intSolPhysData[EulerTerm::VY] - uNX2*normal[YY];
    m_ghostSolPhysData[EulerTerm::VZ]  = m_intSolPhysData[EulerTerm::VZ] - uNX2*normal[ZZ];
    m_ghostSolPhysData[EulerTerm::P]   = m_intSolPhysData[EulerTerm::P];
    m_ghostSolPhysData[EulerTerm::H]   = (gammaDivGammaMinus1*m_ghostSolPhysData[EulerTerm::P]
                                            + 0.5*m_ghostSolPhysData[EulerTerm::RHO]*
                                                  m_intSolPhysData[EulerTerm::V]*
                                                  m_intSolPhysData[EulerTerm::V]
                                         )/m_ghostSolPhysData[EulerTerm::RHO];

    // set the ghost state from its physical data
    m_eulerVarSet->computeStateFromPhysicalData(m_ghostSolPhysData,ghostState);
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCCurvedWallEuler3D::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                                                std::vector< std::vector< RealVector* > >& ghostGrads,
                                                const std::vector< RealVector >& normals,
                                                const std::vector< RealVector >& coords)
{
  throw Common::NotImplementedException (FromHere(),"computeGhostGradients for CurvedWallEuler3D not implemented yet");
}

//////////////////////////////////////////////////////////////////////////////

void BCCurvedWallEuler3D::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  BCStateComputer::configure(args);

  // parsing the functions that the user inputed
  m_vFunction.setFunctions(m_functions);
  m_vFunction.setVariables(m_vars);
  try
  {
    m_vFunction.parse();
  }
  catch (Common::ParserException& e)
  {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCCurvedWallEuler3D::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = false;

  // get Euler 3D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler3DVarSet>();
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler3DVarSet in BCCurvedWallEuler3D!");
  }

  // resize the physical data for internal and ghost solution points
  m_eulerVarSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
  m_eulerVarSet->getModel()->resizePhysicalData(m_intSolPhysData  );

  // get datahandle to nodes
  m_nodes  = MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();

  // Get flux points wheight coordinates
  vector< SpectralFVElementData* >& svLocalData = getMethodData().getSVLocalData();
  std::string computerType = getMethodData().getBndFaceTermComputer()->getType();
  if (computerType == "Std")
  {
    m_flxPntWheightCoordsSVFaces = svLocalData[0]->getExtQPntWheightCoords();
  }
  if (computerType == "QuadFree")
  {
    m_flxPntWheightCoordsSVFaces = svLocalData[0]->getFaceFluxPolyNodeWheightCoord();
  }
  cf_assert(m_flxPntWheightCoordsSVFaces.isNotNull());

  if (m_functions.size() > 0)
  {
    computeNormalsInFaceFluxPointsFromUserFunction();
  }
  else
  {
    computeNormalsInFaceFluxPointsFromMeshData();
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCCurvedWallEuler3D::computeNormalsInFaceFluxPointsFromUserFunction()
{
  CFAUTOTRACE;

  // get TRS list
  vector< SafePtr< TopologicalRegionSet > > trsList = MeshDataStack::getActive()->getTrsList();

  // number of boundary TRSs with this BC
  const CFuint nbrBCTRSs = m_trsNames.size();

  // get boundary TRSs
  vector< SafePtr< TopologicalRegionSet > > m_bcTRSs(nbrBCTRSs);
  const CFuint nbTRSs = trsList.size();
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS)
  {
    for (CFuint iBCTRS = 0; iBCTRS < nbrBCTRSs; ++iBCTRS)
    {
      if (m_trsNames[iBCTRS] == trsList[iTRS]->getName())
      {
        if (m_bcTRSs[iBCTRS].isNull())
        {
          m_bcTRSs[iBCTRS] = trsList[iTRS];
        }
        else
        {
          throw BadFormatException (FromHere(),"Two TRSs with the same name found!");
        }
      }
    }
  }

  // loop over boundary condition TRSs
  for (CFuint iBCTRS = 0; iBCTRS < nbrBCTRSs; ++iBCTRS)
  {
    if (m_bcTRSs[iBCTRS].isNull())
    {
      throw BadFormatException (FromHere(),"Not all boundary TRSs found!");
    }

    // get boundary face connectivities
    SafePtr< ConnectivityTable< CFuint > > bndFaceNodesConn = MeshDataStack::getActive()->getConnectivity(m_trsNames[iBCTRS]+"Nodes");
    SafePtr< vector< CFuint > >            bndFaceLocalIDs  = m_bcTRSs[iBCTRS]->getGeoEntsLocalIdx();

    // loop over faces
    const CFuint nbrFaces = bndFaceNodesConn->nbRows();
    for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
    {
      // get boundary face ID
      const CFuint faceID = (*bndFaceLocalIDs)[iFace];

      // get face nodes
      vector< RealVector > nodeCoords(3,RealVector(3));
      nodeCoords[0] = *m_nodes[(*bndFaceNodesConn)(iFace,0)];
      nodeCoords[1] = *m_nodes[(*bndFaceNodesConn)(iFace,1)];
      nodeCoords[2] = *m_nodes[(*bndFaceNodesConn)(iFace,2)];

      // compute normals in flux points
      m_faceFlxPntNormals[faceID] = vector< RealVector >();
      vector< RealVector >& flxPntNormals = m_faceFlxPntNormals[faceID];
      const CFuint nbrFlxPnts = m_flxPntWheightCoordsSVFaces->size();
      flxPntNormals.resize(nbrFlxPnts);
      for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
      {
        // coordinates
        RealVector flxPntCoor = (*m_flxPntWheightCoordsSVFaces)[iFlx][0]*nodeCoords[0] +
                                (*m_flxPntWheightCoordsSVFaces)[iFlx][1]*nodeCoords[1];
                                (*m_flxPntWheightCoordsSVFaces)[iFlx][2]*nodeCoords[2];

        // compute unit normal
        flxPntNormals[iFlx].resize(3);
        m_vFunction.evaluate(flxPntCoor,flxPntNormals[iFlx]);
//         CF_DEBUG_OBJ(flxPntNormals[iFlx]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCCurvedWallEuler3D::computeNormalsInFaceFluxPointsFromMeshData()
{
  CFAUTOTRACE;

  throw Common::NotImplementedException (FromHere(),"computeNormalsInFaceFluxPointsFromMeshData for CurvedWallEuler3D not implemented yet");

  // get TRS list
  vector< SafePtr< TopologicalRegionSet > > trsList = MeshDataStack::getActive()->getTrsList();

  // number of boundary TRSs with this BC
  const CFuint nbrBCTRSs = m_trsNames.size();

  // get boundary TRSs
  vector< SafePtr< TopologicalRegionSet > > m_bcTRSs(nbrBCTRSs);
  const CFuint nbTRSs = trsList.size();
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS)
  {
    for (CFuint iBCTRS = 0; iBCTRS < nbrBCTRSs; ++iBCTRS)
    {
      if (m_trsNames[iBCTRS] == trsList[iTRS]->getName())
      {
        if (m_bcTRSs[iBCTRS].isNull())
        {
          m_bcTRSs[iBCTRS] = trsList[iTRS];
        }
        else
        {
          throw BadFormatException (FromHere(),"Two TRSs with the same name found!");
        }
      }
    }
  }

  // loop over boundary condition TRSs
  for (CFuint iBCTRS = 0; iBCTRS < nbrBCTRSs; ++iBCTRS)
  {
    if (m_bcTRSs[iBCTRS].isNull())
    {
      throw BadFormatException (FromHere(),"Not all boundary TRSs found!");
    }

    // get boundary face connectivities
    SafePtr< ConnectivityTable< CFuint > > bndFaceNodesConn = MeshDataStack::getActive()->getConnectivity(m_trsNames[iBCTRS]+"Nodes");
    SafePtr< vector< CFuint > >            bndFaceLocalIDs  = m_bcTRSs[iBCTRS]->getGeoEntsLocalIdx();

    // loop over faces
    const CFuint nbrFaces = bndFaceNodesConn->nbRows();
    for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
    {
      // get boundary face ID
      const CFuint faceID = (*bndFaceLocalIDs)[iFace];

      // variables for neighbouring nodes
      vector< bool >       nghbrNodeFound (4,false);
      vector< CFuint >     nghbrNodeIDs   (4);
      vector< RealVector > nghbrNodeCoords(4);

      // put current face nodes
      nghbrNodeFound[1] = true;
      nghbrNodeIDs  [1] = (*bndFaceNodesConn)(iFace,0);
      nghbrNodeFound[2] = true;
      nghbrNodeIDs  [2] = (*bndFaceNodesConn)(iFace,1);

      // find neighbouring face nodes
      for (CFuint iFace2 = 0; iFace2 < nbrFaces && !(nghbrNodeFound[0] && nghbrNodeFound[3]); ++iFace2)
      {
        if (iFace2 != iFace)
        {
          // neighbouring face nodes
          vector< CFuint > nghbrFaceNodes(2);
          nghbrFaceNodes[LEFT ] = (*bndFaceNodesConn)(iFace2,LEFT );
          nghbrFaceNodes[RIGHT] = (*bndFaceNodesConn)(iFace2,RIGHT);

          // if this is a neighbouring face, put the other node in the list
          if (nghbrFaceNodes[LEFT ] == nghbrNodeIDs[1] ||
              nghbrFaceNodes[RIGHT] == nghbrNodeIDs[1])
          {
            if (nghbrNodeFound[0])
            {
              throw BadFormatException (FromHere(),"Face with two neighbouring faces on the same side...");
            }
            nghbrNodeIDs  [0] = nghbrFaceNodes[LEFT ] == nghbrNodeIDs[1] ?
                                nghbrFaceNodes[RIGHT] :
                                nghbrFaceNodes[LEFT ];
            nghbrNodeFound[0] = true;
          }
          if (nghbrFaceNodes[LEFT ] == nghbrNodeIDs[2] ||
              nghbrFaceNodes[RIGHT] == nghbrNodeIDs[2])
          {
            if (nghbrNodeFound[3])
            {
              throw BadFormatException (FromHere(),"Face with two neighbouring faces on the same side...");
            }
            nghbrNodeIDs  [3] = nghbrFaceNodes[LEFT ] == nghbrNodeIDs[2] ?
                                nghbrFaceNodes[RIGHT] :
                                nghbrFaceNodes[LEFT ];
            nghbrNodeFound[3] = true;
          }
        } // if
      } // for

      // set node coordinates
      for (CFuint iNode = 0; iNode < 4; ++iNode)
      {
//         CF_DEBUG_OBJ(nghbrNodeFound[iNode]);
        if (nghbrNodeFound[iNode])
        {
          nghbrNodeCoords[iNode].resize(2);
          nghbrNodeCoords[iNode] = *m_nodes[nghbrNodeIDs[iNode]];
        }
      }

      // compute coordinates for missing nodes
      if (!nghbrNodeFound[0] && !nghbrNodeFound[3])
      // this is a case where there is only one face in the whole TRS
      {
        // linearly extrapolate the coordinates of the neighbouring nodes
        nghbrNodeCoords[0].resize(2);
        nghbrNodeCoords[0] = 2.0*nghbrNodeCoords[1] - nghbrNodeCoords[2];
        nghbrNodeFound [0] = true;
        nghbrNodeCoords[3].resize(2);
        nghbrNodeCoords[3] = 2.0*nghbrNodeCoords[2] - nghbrNodeCoords[1];
        nghbrNodeFound [3] = true;
      }
      else if (!nghbrNodeFound[0])
      {
        // search neighbouring node on the other side and reconstruct coordinates
        // of missing node using a cubic polynomial
        // not implemented, might cause oscillations in the reconstructed geometry
        // coordinates are reconstructed using quadratic polynomial below

        // if no neighbouring node is found (TRS with only two faces)
        // quadratically extrapolate the node coordinates
        if (!nghbrNodeFound[0])
        {
          nghbrNodeCoords[0].resize(2);
          nghbrNodeCoords[0] = 3.0*nghbrNodeCoords[1] - 3.0*nghbrNodeCoords[2] + nghbrNodeCoords[3];
          nghbrNodeFound [0] = true;
        }
      }
      else if (!nghbrNodeFound[3])
      {
        // search neighbouring node on the other side and reconstruct coordinates
        // of missing node using a cubic polynomial
        // not implemented, might cause oscillations in the reconstructed geometry
        // coordinates are reconstructed using quadratic polynomial below

        // if no neighbouring node is found (TRS with only two faces)
        // quadratically extrapolate the node coordinates
        if (!nghbrNodeFound[3])
        {
          nghbrNodeCoords[3].resize(2);
          nghbrNodeCoords[3] = 3.0*nghbrNodeCoords[2] - 3.0*nghbrNodeCoords[1] + nghbrNodeCoords[0];
          nghbrNodeFound [3] = true;
        }
      }

//       CF_DEBUG_OBJ(nghbrNodeCoords[1]);
//       CF_DEBUG_OBJ(nghbrNodeCoords[2]);

      // compute normals in flux points
      m_faceFlxPntNormals[faceID] = vector< RealVector >();
      vector< RealVector >& flxPntNormals = m_faceFlxPntNormals[faceID];
      const CFuint nbrFlxPnts = m_flxPntWheightCoordsSVFaces->size();
      flxPntNormals.resize(nbrFlxPnts);
      RealVector shapeFuncDerivs(4);
      for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
      {
        // flux point local coordinates
        const CFreal ksi = (*m_flxPntWheightCoordsSVFaces)[iFlx][1];
        const CFreal eta = (*m_flxPntWheightCoordsSVFaces)[iFlx][2];

        // get derivative of shape functions at this point
        setShapeFunctionDerivatives(ksi,eta,shapeFuncDerivs);

        // get derivative of absolute coordinates at this point
        const RealVector drdksi = shapeFuncDerivs[0]*nghbrNodeCoords[0] +
                                  shapeFuncDerivs[1]*nghbrNodeCoords[1] +
                                  shapeFuncDerivs[2]*nghbrNodeCoords[2] +
                                  shapeFuncDerivs[3]*nghbrNodeCoords[3];

        // compute unit normal
        flxPntNormals[iFlx].resize(2);
        const CFreal factor = 1.0/sqrt(drdksi[XX]*drdksi[XX] + drdksi[YY]*drdksi[YY]);
        flxPntNormals[iFlx][XX] = +drdksi[YY]*factor;
        flxPntNormals[iFlx][YY] = -drdksi[XX]*factor;

/*        // coordinates
        RealVector flxPntCoor = (*m_flxPntWheightCoordsSVFaces)[iFlx][0]*nghbrNodeCoords[1] +
                                (*m_flxPntWheightCoordsSVFaces)[iFlx][1]*nghbrNodeCoords[2];

        // compute unit normal
        flxPntNormals[iFlx].resize(2);
        const CFreal factor = 1.0/sqrt(flxPntCoor[XX]*flxPntCoor[XX] + flxPntCoor[YY]*flxPntCoor[YY]);
        flxPntNormals[iFlx][XX] = flxPntCoor[XX]*factor;
        flxPntNormals[iFlx][YY] = flxPntCoor[YY]*factor;*/
//         CF_DEBUG_OBJ(flxPntNormals[iFlx]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCCurvedWallEuler3D::setShapeFunctionDerivatives(const CFreal ksi, const CFreal eta,
                                                      RealVector& shapeFuncDerivs)
{
  cf_assert(shapeFuncDerivs.size() == 6);

/*  shapeFuncDerivs[0] = - 1.0/3.0 +     ksi - 0.5*ksi*ksi;
  shapeFuncDerivs[1] = - 0.5     - 2.0*ksi + 1.5*ksi*ksi;
  shapeFuncDerivs[2] = + 1.0     +     ksi - 1.5*ksi*ksi;
  shapeFuncDerivs[3] = - 1.0/6.0           + 0.5*ksi*ksi;*/
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFV

}  // namespace COOLFluiD
