#include "Environment/DirPaths.hh"

#include "Framework/BadFormatException.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "NavierStokes/Euler2DVarSet.hh"
#include "NavierStokes/EulerTerm.hh"

#include "SpectralFD/SpectralFDElementData.hh"

#include "SpectralFDNavierStokes/BCCurvedWallEuler2D.hh"
#include "SpectralFDNavierStokes/SpectralFDNavierStokes.hh"

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

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

Framework::MethodStrategyProvider< BCCurvedWallEuler2D,
                                   SpectralFDMethodData,
                                   BCStateComputer,
                                   SpectralFDNavierStokesModule
                                 > BCCurvedWallEuler2DProvider("CurvedWallEuler2D");

//////////////////////////////////////////////////////////////////////////////

void BCCurvedWallEuler2D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("NormalsDef","Definition of the Functions that describe the normals.");
  options.addConfigOption< bool >("UseDomainModel","Boolean telling whether to use the domain model for the computation of the normals.");
}

//////////////////////////////////////////////////////////////////////////////

BCCurvedWallEuler2D::BCCurvedWallEuler2D(const std::string& name) :
  BCStateComputer(name),
  m_nodes(CFNULL),
  m_eulerVarSet(CFNULL),
  m_ghostSolPhysData(),
  m_intSolPhysData(),
  m_flxPntFaceLocalCoords(CFNULL),
  m_faceFlxPntNormals(),
  m_useDomainModel()
{
  addConfigOptionsTo(this);

  m_functions = vector<std::string>();
  setParameter("NormalsDef",&m_functions);

  m_vars = vector<std::string>();
  setParameter("Vars",&m_vars);

  m_useDomainModel = false;
  setParameter("UseDomainModel",&m_useDomainModel);
}

//////////////////////////////////////////////////////////////////////////////

BCCurvedWallEuler2D::~BCCurvedWallEuler2D()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCCurvedWallEuler2D::computeGhostStates(const vector< State* >& intStates,
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
  const CFreal gammaDivGammaMinus1 = gamma/(gamma -1.0);

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

    cf_assert(intState.size() == 4);
    cf_assert(ghostState.size() == 4);

    // set the physical data starting from the inner state
    m_eulerVarSet->computePhysicalData(intState,m_intSolPhysData);

    // compute normal velocity component
    const CFreal uNX2 = 2.0*(m_intSolPhysData[EulerTerm::VX]*normal[XX] +
        m_intSolPhysData[EulerTerm::VY]*normal[YY]);

    // set the physical data for the ghost state
    m_ghostSolPhysData[EulerTerm::RHO] = m_intSolPhysData[EulerTerm::RHO];
    m_ghostSolPhysData[EulerTerm::VX]  = m_intSolPhysData[EulerTerm::VX] - uNX2*normal[XX];
    m_ghostSolPhysData[EulerTerm::VY]  = m_intSolPhysData[EulerTerm::VY] - uNX2*normal[YY];
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

void BCCurvedWallEuler2D::computeGhostGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                                                std::vector< std::vector< RealVector* > >& ghostGrads,
                                                const std::vector< RealVector >& normals,
                                                const std::vector< RealVector >& coords)
{
  // number of state gradients
  const CFuint nbrStateGrads = intGrads.size();
  cf_assert(nbrStateGrads == ghostGrads.size());
  cf_assert(nbrStateGrads == normals.size());

  // number of gradient variables
  cf_assert(nbrStateGrads > 0);

  // get the real normals for the current face
  vector< RealVector >& flxPntNormals = m_faceFlxPntNormals[m_faceID];

  // set the ghost gradients
  for (CFuint iState = 0; iState < nbrStateGrads; ++iState)
  {
    cf_assert(intGrads[iState].size() == 4);

    // normal
    const RealVector& normal = flxPntNormals[iState];

    // tangential unit vector
    RealVector tangent(2);
    tangent[XX] = -normal[YY];
    tangent[YY] =  normal[XX];

    // pressure
    RealVector& presGradI = *intGrads  [iState][0];
    RealVector& presGradG = *ghostGrads[iState][0];
    const CFreal nPresGrad = (presGradI[XX]*normal[XX] + presGradI[YY]*normal[YY]);
    presGradG = presGradI - 2.0*nPresGrad*normal;

    // temperature
    RealVector& tempGradI = *intGrads  [iState][3];
    RealVector& tempGradG = *ghostGrads[iState][3];
    const CFreal nTempGrad = (tempGradI[XX]*normal[XX] + tempGradI[YY]*normal[YY]);
    tempGradG = tempGradI - 2.0*nTempGrad*normal;

    // velocity
    RealVector& uGradI = *intGrads  [iState][1];
    RealVector& uGradG = *ghostGrads[iState][1];
    RealVector& vGradI = *intGrads  [iState][2];
    RealVector& vGradG = *ghostGrads[iState][2];

    // internal normal and tangential component
    const RealVector uNGradI = uGradI*normal [XX] + vGradI*normal [YY];
    const RealVector uTGradI = uGradI*tangent[XX] + vGradI*tangent[YY];

    // ghost normal and tangential component
    const RealVector uNGradG = uNGradI;
    const CFreal nGradUT = uTGradI[XX]*normal[XX] + uTGradI[YY]*normal[YY];
    const RealVector uTGradG = uTGradI - 2.0*nGradUT*normal;

    // project onto x- and y-axis
    uGradG = uNGradG*normal[XX] + uTGradG*tangent[XX];
    vGradG = uNGradG*normal[YY] + uTGradG*tangent[YY];
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCCurvedWallEuler2D::configure ( Config::ConfigArgs& args )
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

void BCCurvedWallEuler2D::setup()
{
  CFAUTOTRACE;

  // setup of the parent class
  BCStateComputer::setup();

  // no flux point coordinates required
  m_needsSpatCoord = false;

  // get Euler 2D varset
  m_eulerVarSet = getMethodData().getUpdateVar().d_castTo<Euler2DVarSet>();
  if (m_eulerVarSet.isNull())
  {
    throw Common::ShouldNotBeHereException (FromHere(),"Update variable set is not Euler2DVarSet in BCCurvedWallEuler2D!");
  }

  // resize the physical data for internal and ghost solution points
  m_eulerVarSet->getModel()->resizePhysicalData(m_ghostSolPhysData);
  m_eulerVarSet->getModel()->resizePhysicalData(m_intSolPhysData  );

  // get datahandle to nodes
  m_nodes  = MeshDataStack::getActive()->getNodeDataSocketSink().getDataHandle();

  // Get flux points wheight coordinates
  vector< SpectralFDElementData* >& sdLocalData = getMethodData().getSDLocalData();
  m_flxPntFaceLocalCoords = sdLocalData[0]->getFaceFlxPntsFaceLocalCoords();
  cf_assert(m_flxPntFaceLocalCoords.isNotNull());

  // compute the normals in the face flux points
  if (m_functions.size() > 0)
  {
    computeNormalsInFaceFluxPointsFromUserFunction();
  }
  else if (m_useDomainModel)
  {
    computeNormalsInFaceFluxPointsFromDomainModel();
  }
  else
  {
    computeNormalsInFaceFluxPointsFromMeshData();
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCCurvedWallEuler2D::computeNormalsInFaceFluxPointsFromUserFunction()
{
  CFAUTOTRACE;

  // get Inner Cells TRS
  SafePtr<TopologicalRegionSet> cellTRS = MeshDataStack::getActive()->getTrs("InnerCells");

  // get face builder and set cells trs
  SafePtr< GeometricEntityPool<FaceToCellGEBuilder> > geoBuilder = getMethodData().getFaceBuilder();
  FaceToCellGEBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.cellsTRS = cellTRS;

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

    // set face trs
    geoData.facesTRS = m_bcTRSs[iBCTRS];
    geoData.isBoundary = true;

    // loop over faces
    const CFuint nbrFaces = m_bcTRSs[iBCTRS]->getLocalNbGeoEnts();
    for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
    {
      // build the GeometricEntity
      geoData.idx = iFace;
      GeometricEntity& face = *geoBuilder->buildGE();

      // get boundary face ID
      const CFuint faceID = face.getID();

      // compute normals in flux points
      m_faceFlxPntNormals[faceID] = vector< RealVector >();
      vector< RealVector >& flxPntNormals = m_faceFlxPntNormals[faceID];
      const CFuint nbrFlxPnts = m_flxPntFaceLocalCoords->size();
      flxPntNormals.resize(nbrFlxPnts);
      for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
      {
        // coordinates
        const RealVector flxPntCoor
            = face.computeCoordFromMappedCoord((*m_flxPntFaceLocalCoords)[iFlx]);

        // compute unit normal
        flxPntNormals[iFlx].resize(2);
        m_vFunction.evaluate(flxPntCoor,flxPntNormals[iFlx]);
      }

      // release the GeometricEntity
      geoBuilder->releaseGE();
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCCurvedWallEuler2D::computeNormalsInFaceFluxPointsFromMeshData()
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
      const CFuint nbrFlxPnts = m_flxPntFaceLocalCoords->size();
      flxPntNormals.resize(nbrFlxPnts);
      RealVector shapeFuncDerivs(4);
      for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
      {
        // flux point local coordinate
        const CFreal ksi = (*m_flxPntFaceLocalCoords)[iFlx][KSI];

        // get derivative of shape functions at this point
        setShapeFunctionDerivatives(ksi,shapeFuncDerivs);

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

//         CF_DEBUG_OBJ(flxPntNormals[iFlx]);
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCCurvedWallEuler2D::computeNormalsInFaceFluxPointsFromDomainModel()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCCurvedWallEuler2D::setShapeFunctionDerivatives(const CFreal ksi, RealVector& shapeFuncDerivs)
{
  cf_assert(shapeFuncDerivs.size() == 4);

  shapeFuncDerivs[0] = - 1.0/3.0 +     ksi - 0.5*ksi*ksi;
  shapeFuncDerivs[1] = - 0.5     - 2.0*ksi + 1.5*ksi*ksi;
  shapeFuncDerivs[2] = + 1.0     +     ksi - 1.5*ksi*ksi;
  shapeFuncDerivs[3] = - 1.0/6.0           + 0.5*ksi*ksi;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD
