#include "Framework/MethodCommandProvider.hh"
#include "Framework/CFSide.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
// #include "Environment/SingleBehaviorFactory.hh"

#include "DiscontGalerkin/StdStabilization.hh"
#include "DiscontGalerkin/DiscontGalerkin.hh"
#include "Framework/FaceToCellGEBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StdStabilization,DiscontGalerkinSolverData,DiscontGalerkinModule >
  StdStabilizationProvider("StdStabilization");

//////////////////////////////////////////////////////////////////////////////

StdStabilization::StdStabilization(const std::string& name)
  : StdBaseSolve(name),
    m_cells(CFNULL),
    m_mapElemData(),
    socket_integrationIndex("integrationIndex"),
    socket_normals("normals"),
    m_gOnMinusOne(),
    m_diameter(),
    m_gRhoJump()
{
}

//////////////////////////////////////////////////////////////////////////////

StdStabilization::~StdStabilization()
{
}

//////////////////////////////////////////////////////////////////////////////

void StdStabilization::setup()
{
  CFAUTOTRACE;
  StdBaseSolve::setup();
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  // set a pointer to the cells and faces
  m_cells.reset(MeshDataStack::getActive()->getTrs("InnerCells"));

  //get types of elements in triangulation
  SafePtr<vector<ElementTypeData> > elementType = MeshDataStack::getActive()->getElementTypeData();


  const CFuint nbElemTypes = elementType->size();

  // loop over types since it can happen to deal with an hybrid mesh
  for (CFuint iType = 0; iType < nbElemTypes; ++iType) {

     const CFuint nbStatesInType = (*elementType)[iType].getNbStates();

     BlockAccumulator* ptr = getMethodData().getLinearSystemSolver()[0]->
       createBlockAccumulator(nbStatesInType,nbStatesInType,nbEqs);
     RealVector* vec = new RealVector(nbStatesInType*nbEqs);
     RealVector* vec2 = new RealVector(nbEqs);
     RealMatrix* mat = new RealMatrix(nbStatesInType*nbEqs,nbStatesInType*nbEqs);
     vector<RealVector>* residual = new vector<RealVector>(nbStatesInType,*vec2);

     DGElemTypeData elemTypeData(ptr,mat,vec,residual);

     m_mapElemData.insert(nbStatesInType, elemTypeData);
  }

  // loop over types since it can happen to deal with an hybrid mesh
  for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
    for (CFuint jType = iType; jType < nbElemTypes; ++jType) {

      const CFuint nbStatesInType = (*elementType)[iType].getNbStates()+(*elementType)[jType].getNbStates();

      BlockAccumulator* ptr = getMethodData().getLinearSystemSolver()[0]->
       createBlockAccumulator(nbStatesInType,nbStatesInType,nbEqs);
      RealVector* vec = new RealVector(nbStatesInType*nbEqs);
      RealVector* vec2 = new RealVector(nbEqs);
      RealMatrix* mat = new RealMatrix(nbStatesInType*nbEqs,nbStatesInType*nbEqs);
      vector<RealVector>* residual = new vector<RealVector>(nbStatesInType,*vec2);

      DGElemTypeData elemTypeData(ptr,mat,vec,residual);

      m_mapElemData.insert(nbStatesInType, elemTypeData);
    }
  }

  m_mapElemData.sortKeys();

  //make help vectors of diameters and gOnMinusOne
  SafePtr<TopologicalRegionSet> trs = m_cells;

  //get number of inner cells in region inner cells
  const CFuint nbGeos = trs->getLocalNbGeoEnts();

  m_diameter.resize(nbGeos);
  m_gOnMinusOne.resize(nbGeos);
  m_gRhoJump.resize(nbGeos);
  //get geobuilder to build cells with needed properties (connection, ..)
  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();


  //get structure of data from geobuilder (cells with properties)
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  //set that we use only data of inner cells trs
  geoData.trs = trs;

  //loop over all inner cells
  for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt) 
  {
    //set index of cell
    geoData.idx = iGeoEnt;
    //geo builder make cell 

    GeometricEntity& cell = *geoBuilder->buildGE();

    //read vector of nodes on element (cell)
    std::vector<Node*>&  nodes  = *cell.getNodes();

    m_diameter[iGeoEnt] = 0.0;
    //only for TRIAG and TETRA
    if (nbDim == 2)
    {
      for(CFuint i=0;i<3;i++)
      {
        for(CFuint j=i+1;j<3;j++)
        {
          CFreal diam=sqrt(pow(((*nodes[i])[YY]-(*nodes[j])[YY]),2) + pow(((*nodes[i])[XX]-(*nodes[j])[XX]),2));
          if (diam > m_diameter[iGeoEnt]) m_diameter[iGeoEnt]=diam;
        }
      }
      m_gOnMinusOne[iGeoEnt]=0.0;
      for(CFuint i=0;i<3;i++)
      {
        m_gOnMinusOne[iGeoEnt] += sqrt(pow((*nodes[(i+1)%3])[0] - (*nodes[i])[0],2) + pow((*nodes[(i+1)%3])[1] - (*nodes[i])[1],2));
      }
      //computation of the Jacobi determinant of mapping from refference element to cell
      m_gOnMinusOne[iGeoEnt] *= pow((abs(cell.computeVolume())*2.0),3.0/4.0);
    }
    else
    {
      for(CFuint i=0;i<4;i++)
      {
        for(CFuint j=i;j<4;j++)
        {
          CFreal diam=sqrt(pow(((*nodes[i])[ZZ]-(*nodes[j])[ZZ]),2) + pow(((*nodes[i])[YY]-(*nodes[j])[YY]),2) + pow(((*nodes[i])[XX]-(*nodes[j])[XX]),2));
          if (diam > m_diameter[iGeoEnt]) m_diameter[iGeoEnt]=diam;
        }
      }
      m_gOnMinusOne[iGeoEnt]=0.0;
      RealVector vect1(3);
      RealVector vect2(3);

      for(CFuint i=0;i<3;i++)
      {
        for(CFuint j=0;j<3;j++)
        {
          vect1[j]= (*nodes[(i)+1])[j]-(*nodes[0])[j];
          vect2[j]= (*nodes[(i+1)%3 + 1])[j]-(*nodes[0])[j];
        }
        m_gOnMinusOne[iGeoEnt] += sqrt(pow(vect1[0]*vect2[1]-vect1[1]*vect2[0],2) + pow(vect1[1]*vect2[2]-vect1[2]*vect2[1],2) + pow(vect1[0]*vect2[2]-vect1[2]*vect2[0],2));
      }
      for(CFuint j=0;j<3;j++)
      {
        vect1[j]= (*nodes[2])[j]-(*nodes[1])[j];
        vect2[j]= (*nodes[3])[j]-(*nodes[1])[j];
      }
      m_gOnMinusOne[iGeoEnt] += sqrt(pow(vect1[0]*vect2[1]-vect1[1]*vect2[0],2) + pow(vect1[1]*vect2[2]-vect1[2]*vect2[1],2) + pow(vect1[0]*vect2[2]-vect1[2]*vect2[0],2));
      //computation of the Jacobi determinant of mapping from refference element to cell
      m_gOnMinusOne[iGeoEnt] *= pow((abs(cell.computeVolume())*6.0),1.0/6.0);
    }
    geoBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdStabilization::unsetup()
{
  CFAUTOTRACE;
  StdBaseSolve::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void StdStabilization::execute()
{
  CFAUTOTRACE;
  CFout << "StdStabilization is applied " << CFendl;
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  //computation of g_Ki
  CFreal rho;

  DataHandle< std::vector< CFuint > >
    integrationIndex = socket_integrationIndex.getDataHandle();
  DataHandle< RealVector > normals = socket_normals.getDataHandle();
  Common::SafePtr<GeometricEntityPool<Framework::FaceToCellGEBuilder> > geoBuilderFace = getMethodData().getFaceBuilder();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get InnerFaces TopologicalRegionSet
  SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");

  //get number of inner cells in region inner cells
  const CFuint nbCells = cells->getLocalNbGeoEnts();

  for(CFuint i=0;i < nbCells; i++)
  {
    m_gRhoJump[i]=0.0;
  }

  // get the geodata of the face builder and set the TRSs
  Framework::FaceToCellGEBuilder::GeoData& geoDataFace = geoBuilderFace->getDataGE();
  geoDataFace.cellsTRS = cells;
  geoDataFace.facesTRS = faces;
  geoDataFace.isBoundary = false;

  const CFuint nbFaces = faces->getLocalNbGeoEnts();
  //loop over all inner faces

  for (CFuint iFace = 0; iFace < nbFaces; ++iFace)
  {
    //set index of face
    geoDataFace.idx = iFace;
    //geo builder make face
    GeometricEntity& face = *geoBuilderFace->buildGE();

    //nodes of the actual face
    std::vector<Node*>&  nodes  = *face.getNodes();
    CFreal detJacobi;
    if (nbDim == 2)
    {
      RealVector hlp = *nodes[1] - *nodes[0];
      detJacobi = sqrt(hlp[0]*hlp[0]+hlp[1]*hlp[1]);
    }
    else
    {
      //     detJacobi = face.computeVolume()*2;
      RealVector hlp1 = *nodes[1] - *nodes[0];
      RealVector hlp2 = *nodes[2] - *nodes[0];
      detJacobi = sqrt((hlp1[1]*hlp2[2] - hlp1[2]*hlp2[1])*(hlp1[1]*hlp2[2] - hlp1[2]*hlp2[1]) + (hlp1[2]*hlp2[0] - hlp1[0]*hlp2[2])*(hlp1[2]*hlp2[0] - hlp1[0]*hlp2[2])+(hlp1[0]*hlp2[1] - hlp1[1]*hlp2[0])*(hlp1[0]*hlp2[1] - hlp1[1]*hlp2[0]));//*2.0/2.0;
    }

    //get the left (neighbouring) cell
    GeometricEntity* cellLeft  = face.getNeighborGeo(LEFT);

    //get nodes and states of the left cell
//     std::vector<Node*>&  left_cell_nodes  = *cellLeft->getNodes();
    std::vector<State*>& left_cell_states = *cellLeft->getStates();

    //get the right (neighbouring) cell
    GeometricEntity* cellRight = face.getNeighborGeo(RIGHT);
    //get nodes and states of the right cell
//     std::vector<Node*>&  right_cell_nodes  = *cellRight->getNodes();
    std::vector<State*>& right_cell_states = *cellRight->getStates();

    //get number of states in left and right cell
    const CFuint nbStatesInCellLeft  = left_cell_states.size();
    const CFuint nbStatesInCellRight = right_cell_states.size();

    //*****************************************************************
    //*****************************************************************
    //SET Integrator
    //compute shape function in quadrature points
    const std::vector<RealVector>& leftShapeFunctions =  getMethodData().getContourIntegrator()->getSolutionIntegrator(cellLeft)->computeShapeFunctionsAtQuadraturePoints();

    //compute shape function in quadrature points
    const std::vector<RealVector>& rightShapeFunctions =  getMethodData().getContourIntegrator()->getSolutionIntegrator(cellRight)->computeShapeFunctionsAtQuadraturePoints();


    //numbers of quadrature points
    CFuint m_nbKvadrPoint = getMethodData().getContourIntegrator()->getSolutionIntegrator(cellLeft)->getIntegratorPattern()[0];

    //set weights for element quadrature
    const std::vector<RealVector>& leftWeight = getMethodData().getContourIntegrator()->getSolutionIntegrator(cellLeft)->getCoeff();

//     const std::vector<RealVector>& rightWeight = getMethodData().getContourIntegrator()->getSolutionIntegrator(cellRight)->getCoeff();

    //get coordinates of quadrature points
    const std::vector<RealVector>& leftCoord = getMethodData().getContourIntegrator()->getSolutionIntegrator(cellLeft)->getQuadraturePointsCoordinates();

    //get coordinates of quadrature points
    const std::vector<RealVector>& rightCoord = getMethodData().getContourIntegrator()->getSolutionIntegrator(cellRight)->getQuadraturePointsCoordinates();

    //compute gradient of shape functions in quadrature points
    std::vector<RealMatrix> leftGradient = cellLeft->computeSolutionShapeFunctionGradientsInMappedCoordinates(leftCoord);


    //compute gradient of shape functions in quadrature points
    std::vector<RealMatrix> rightGradient = cellRight->computeSolutionShapeFunctionGradientsInMappedCoordinates(rightCoord);

    //*****************************************************************
    //*****************************************************************

    RealVector normal(nbDim);
    for(CFuint i=0;i<nbDim;i++) normal[i]=normals[iFace][i];
    CFuint m_idxFaceFromLeftCell=integrationIndex[iFace][0];
    CFuint m_idxFaceFromRightCell=integrationIndex[iFace][1];
    CFuint rightHlpIndex=integrationIndex[iFace][2];
    CFuint swifted = integrationIndex[iFace][3];

    CFreal temp_value=0.0;
    for(CFuint kvadrature_point = 0; kvadrature_point < m_nbKvadrPoint; kvadrature_point++ )
    {
      CFuint leftIndex  = m_idxFaceFromLeftCell*m_nbKvadrPoint + kvadrature_point;
      CFuint rightIndex;
      if (nbDim == 2)
      {
        rightIndex = m_idxFaceFromRightCell*m_nbKvadrPoint + m_nbKvadrPoint - 1 - kvadrature_point;
      }
      else
      {
        if (kvadrature_point!=0)
        {
          if (swifted == 1)
          {
            rightIndex = m_idxFaceFromRightCell*m_nbKvadrPoint + 1 + ((kvadrature_point-1)/3)*3 + (2-(7+rightHlpIndex-kvadrature_point +1)%3);
          }
          else
          {
           rightIndex = m_idxFaceFromRightCell*m_nbKvadrPoint + 1 + ((kvadrature_point-1)/3)*3 + (2-(kvadrature_point+rightHlpIndex)%3);
          }
        }
        else
        {
          rightIndex = m_idxFaceFromRightCell*m_nbKvadrPoint + kvadrature_point;
        }
      }
      RealVector hlpVector=cellLeft->computeCoordFromMappedCoord(leftCoord[leftIndex]) - cellRight->computeCoordFromMappedCoord(rightCoord[rightIndex]);

      assert((hlpVector.norm1()/detJacobi < 0.0001));

      //computation of jump of rho
      rho =0.0;
      for (CFuint iState = 0; iState < nbStatesInCellLeft; ++iState) //loop over states in cell
      {
        RealVector &states = *left_cell_states[iState]->getData();
        rho += leftShapeFunctions[leftIndex][iState]*states[0];
      }
      for (CFuint iState = 0; iState < nbStatesInCellRight; ++iState) //loop over states in cell 
      {
        RealVector &states = *right_cell_states[iState]->getData();
        rho -= rightShapeFunctions[rightIndex][iState]*states[0];
      }
      rho *= rho;

       temp_value += rho*leftWeight[0][kvadrature_point];

    }
    temp_value *=detJacobi;
    m_gRhoJump[cellLeft->getID()] += temp_value;
    m_gRhoJump[cellRight->getID()] += temp_value;

    // release the face
    geoBuilderFace->releaseGE();
  }


//computation of stabilization term on elements - volume integrators
CFout << "to cells" << CFendl;
//get matrix of linear systemsolver
  SafePtr<LinearSystemSolver> lss = getMethodData().getLinearSystemSolver()[0];
  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();

  //get geobuilder to build cells with needed properties (connection, ..)
  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilderCell = getMethodData().getStdTrsGeoBuilder();


  //get structure of data from geobuilder (cells with properties)
  StdTrsGeoBuilder::GeoData& geoDataCell = geoBuilderCell->getDataGE();
  //set that we use only data of inner cells trs
  geoDataCell.trs =cells;

  //loop over all inner cells
  for(CFuint iCell = 0; iCell < nbCells; ++iCell) {
    CFLogDebugMax("Cell " << iCell << "\n");

    // build the GeometricEntity (cell)
    //set index of cell
    geoDataCell.idx = iCell;
    //geo builder make cell 
    GeometricEntity& cell = *geoBuilderCell->buildGE();

//*****************************************************************
//*****************************************************************
    //SET Integrator
    //compute shape function in quadrature points
//     const std::vector<RealVector>& shapeFunctions =  getMethodData().getVolumeIntegrator()->getSolutionIntegrator(&cell)->computeShapeFunctionsAtQuadraturePoints();

    //numbers of quadrature points
    CFuint m_nbKvadrPoint = getMethodData().getVolumeIntegrator()->getSolutionIntegrator(&cell)->getIntegratorPattern()[0];

    //set weights for element quadrature
    const std::valarray<CFreal>& weight = getMethodData().getVolumeIntegrator()->getSolutionIntegrator(&cell)->getCoeff();

    //get coordinates of quadrature points
    const std::vector<RealVector>& coord = getMethodData().getVolumeIntegrator()->getSolutionIntegrator(&cell)->getQuadraturePointsCoordinates();

    //compute gradient of shape functions in quadrature points
    std::vector<RealMatrix> gradient = cell.computeSolutionShapeFunctionGradientsInMappedCoordinates(coord);

//*****************************************************************
//*****************************************************************

    //read vector of states on element (cell)
    std::vector<State*>& cellStates = *(cell.getStates());
    //read number of states on element (cell)
    const CFuint nbStatesInCell = cellStates.size();
    // CFout << nbStatesInCell << "  " << "\n" << CFendl;
    //read vector of nodes on element (cell)

    /// @todo this must be improved. Finding in a map
    /// for every cell is efficiently speaking not acceptable.
    //make block acumulator with size of sum of states in cell
    DGElemTypeData elemData = m_mapElemData.find(nbStatesInCell);
    BlockAccumulator& acc = *elemData.first;
    RealMatrix& elemMat = *elemData.second;

    //set matrix in blockaccumulator to 0
    elemMat=0.0;
    acc.setValuesM(elemMat);

    //computation of the Jacobi determinant of mapping from refference element to cell
    ///this must be generalized
    if (nbDim == 2)
    {
      detJacobi = abs(cell.computeVolume())*2.0;
    }
    else
    {
      detJacobi = abs(cell.computeVolume())*6.0;
    }

    // set the IDs on the blockaccumulator (we use setRowColIndex() )
    //connection between local and global state ID
    for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
      const CFuint stateID = cellStates[iState]->getLocalID();
      acc.setRowColIndex(iState, stateID);
    }
    cf_assert(cell.getID()==iCell);
    //loop over kvadrature point on the cell
    for(CFuint kvadrature_point = 0; kvadrature_point < m_nbKvadrPoint; kvadrature_point++ )
    {
      //set elemMat to 0 if isn't
      if (kvadrature_point!=0) elemMat=0.0;
      //compute inner face term 
      //loop over test function
      for(CFuint row = 0; row < nbStatesInCell; row++ )
      {
        //loop over base function of solution
        for(CFuint col = 0; col < nbStatesInCell; col++ )
        {
          CFreal temp_value=0.0;
          for(CFuint s = 0; s < nbDim; s++ )
          {
            temp_value+=gradient[kvadrature_point](row,s)*gradient[kvadrature_point](col,s);
          }
          for(CFuint i = 0; i < nbEqs; i++ )
          {
            elemMat(row*nbEqs + i, col*nbEqs + i) += temp_value;
          }
        }
      }

      //finaly multiply by quadrature weight
      elemMat*=1.0/m_gOnMinusOne[iCell]*m_gRhoJump[iCell]*weight[kvadrature_point]*detJacobi;
      // add local matrix to matrix of linear solver using block accumulator

      acc.addValuesM(elemMat);
    }
// CFDEBUGOBJ(elemVec);
// CFDEBUGOBJ(elemMat);

    // add the values in the jacobian matrix
    jacobMatrix->addValues(acc);
    //release the GeometricEntity
    geoBuilderCell->releaseGE();
  }


  //computation of stabilization term on faces - contour integrators
  CFout << ", to faces \n" << CFendl;
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace)
  {

    CFLogDebugMax("Face " << iFace << "\n");
    //set index of face
    geoDataFace.idx = iFace;
    //geo builder make face
    GeometricEntity& face = *geoBuilderFace->buildGE();

    //nodes of the actual face
    std::vector<Node*>&  nodes  = *face.getNodes();
    CFreal detJacobi;
    if (nbDim == 2)
    {
      RealVector hlp = *nodes[1] - *nodes[0];
      detJacobi = sqrt(hlp[0]*hlp[0]+hlp[1]*hlp[1]);
    }
    else
    {
      //     detJacobi = face.computeVolume()*2;
      RealVector hlp1 = *nodes[1] - *nodes[0];
      RealVector hlp2 = *nodes[2] - *nodes[0];
      detJacobi = sqrt((hlp1[1]*hlp2[2] - hlp1[2]*hlp2[1])*(hlp1[1]*hlp2[2] - hlp1[2]*hlp2[1]) + (hlp1[2]*hlp2[0] - hlp1[0]*hlp2[2])*(hlp1[2]*hlp2[0] - hlp1[0]*hlp2[2])+(hlp1[0]*hlp2[1] - hlp1[1]*hlp2[0])*(hlp1[0]*hlp2[1] - hlp1[1]*hlp2[0]));//*2.0/2.0;
    }

    //get the left (neighbouring) cell
    GeometricEntity* cellLeft  = face.getNeighborGeo(LEFT);
    CFuint leftID = cellLeft->getID();
    //get nodes and states of the left cell
//     std::vector<Node*>&  left_cell_nodes  = *cellLeft->getNodes();
    std::vector<State*>& left_cell_states = *cellLeft->getStates();

    //get the right (neighbouring) cell
    GeometricEntity* cellRight = face.getNeighborGeo(RIGHT);
    CFuint rightID = cellRight->getID();
    //get nodes and states of the right cell
//     std::vector<Node*>&  right_cell_nodes  = *cellRight->getNodes();
    std::vector<State*>& right_cell_states = *cellRight->getStates();

    //get number of states in left and right cell
    const CFuint nbStatesInCellLeft  = left_cell_states.size();
    const CFuint nbStatesInCellRight = right_cell_states.size();

//*****************************************************************
//*****************************************************************
    //SET Integrator
    //compute shape function in quadrature points
    const std::vector<RealVector>& leftShapeFunctions =  getMethodData().getContourIntegrator()->getSolutionIntegrator(cellLeft)->computeShapeFunctionsAtQuadraturePoints();

    //compute shape function in quadrature points
    const std::vector<RealVector>& rightShapeFunctions =  getMethodData().getContourIntegrator()->getSolutionIntegrator(cellRight)->computeShapeFunctionsAtQuadraturePoints();


    //numbers of quadrature points
    CFuint m_nbKvadrPoint = getMethodData().getContourIntegrator()->getSolutionIntegrator(cellLeft)->getIntegratorPattern()[0];

    //set weights for element quadrature
    const std::vector<RealVector>& leftWeight = getMethodData().getContourIntegrator()->getSolutionIntegrator(cellLeft)->getCoeff();

//     const std::vector<RealVector>& rightWeight = getMethodData().getContourIntegrator()->getSolutionIntegrator(cellRight)->getCoeff();

    //get coordinates of quadrature points
    const std::vector<RealVector>& leftCoord = getMethodData().getContourIntegrator()->getSolutionIntegrator(cellLeft)->getQuadraturePointsCoordinates();
  
    //get coordinates of quadrature points
    const std::vector<RealVector>& rightCoord = getMethodData().getContourIntegrator()->getSolutionIntegrator(cellRight)->getQuadraturePointsCoordinates();

    //compute gradient of shape functions in quadrature points
    std::vector<RealMatrix> leftGradient = cellLeft->computeSolutionShapeFunctionGradientsInMappedCoordinates(leftCoord);


    //compute gradient of shape functions in quadrature points
    std::vector<RealMatrix> rightGradient = cellRight->computeSolutionShapeFunctionGradientsInMappedCoordinates(rightCoord);

//*****************************************************************
//*****************************************************************

    //get number of equations
    const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  
    /// @todo this must be improved. Finding in a map
    /// for every cell is efficiently speaking not acceptable.
    //make block acumulator with size of sum states from left and right cell
    DGElemTypeData elemData = m_mapElemData.find(nbStatesInCellLeft + nbStatesInCellRight);
    BlockAccumulator& acc = *elemData.first;
    RealMatrix& elemMat = *elemData.second;

    //set matrix in blockaccumulator to 0
    elemMat=0.0;
    acc.setValuesM(elemMat);

    // set the IDs on the blockaccumulator (use setRowColIndex() )
    for (CFuint iState = 0; iState < nbStatesInCellLeft; ++iState) {
      const CFuint stateID = left_cell_states[iState]->getLocalID();
      acc.setRowColIndex(iState, stateID);
    }
    for (CFuint iState = 0; iState < nbStatesInCellRight; ++iState) {
      const CFuint stateID = right_cell_states[iState]->getLocalID();
      acc.setRowColIndex(iState + nbStatesInCellLeft, stateID);
    }

    RealVector normal(nbDim);
    for(CFuint i=0;i<nbDim;i++) normal[i]=normals[iFace][i];
    CFuint m_idxFaceFromLeftCell=integrationIndex[iFace][0];
    CFuint m_idxFaceFromRightCell=integrationIndex[iFace][1];
    CFuint rightHlpIndex=integrationIndex[iFace][2];
    CFuint swifted = integrationIndex[iFace][3];

    //loop over kvadrature point on the face
    for(CFuint kvadrature_point = 0; kvadrature_point < m_nbKvadrPoint; kvadrature_point++ )
    {
      CFuint leftIndex  = m_idxFaceFromLeftCell*m_nbKvadrPoint + kvadrature_point;
      CFuint rightIndex;
      if (nbDim == 2)
      {
        rightIndex = m_idxFaceFromRightCell*m_nbKvadrPoint + m_nbKvadrPoint - 1 - kvadrature_point;
      }
      else
      {
        if (kvadrature_point!=0)
        {
          if (swifted == 1)
          {
            rightIndex = m_idxFaceFromRightCell*m_nbKvadrPoint + 1 + ((kvadrature_point-1)/3)*3 + (2-(7+rightHlpIndex-kvadrature_point +1)%3);
          }
          else
          {
           rightIndex = m_idxFaceFromRightCell*m_nbKvadrPoint + 1 + ((kvadrature_point-1)/3)*3 + (2-(kvadrature_point+rightHlpIndex)%3);
          }
        }
        else
        {
          rightIndex = m_idxFaceFromRightCell*m_nbKvadrPoint + kvadrature_point;
        }
      }
      RealVector hlpVector=cellLeft->computeCoordFromMappedCoord(leftCoord[leftIndex]) - cellRight->computeCoordFromMappedCoord(rightCoord[rightIndex]);
      assert((hlpVector.norm1()/detJacobi < 0.0001));

      elemMat=0.0;
      //compute inner face term 
      //loop over test function
      for(CFuint row = 0; row < nbStatesInCellLeft + nbStatesInCellRight; row++ )
      {
        //loop over base function of solution
        for(CFuint col = 0; col < nbStatesInCellLeft + nbStatesInCellRight; col++ )
        {
          if ((col < nbStatesInCellLeft)&&(row < nbStatesInCellLeft))
          {
            //temp_value = multiplication of test function and base function in kvadrature point
            CFreal temp_value=leftShapeFunctions[leftIndex][col]*leftShapeFunctions[leftIndex][row];
            for(CFuint i = 0; i < nbEqs; i++ )
            {
              elemMat(row*nbEqs + i, col*nbEqs + i)+=temp_value;
            }
          }
          //test function is from left cell and base functin from right cell
          else if ((col >= nbStatesInCellLeft)&&(row < nbStatesInCellLeft))
          {
            CFreal temp_value=rightShapeFunctions[rightIndex][col-nbStatesInCellLeft]*leftShapeFunctions[leftIndex][row];
            for(CFuint i = 0; i < nbEqs; i++ )
            {
              elemMat(row*nbEqs + i, col*nbEqs + i)-=temp_value;
            }
           }
          else if ((col < nbStatesInCellLeft)&&(row >= nbStatesInCellLeft))
          {

            CFreal temp_value=rightShapeFunctions[rightIndex][row-nbStatesInCellLeft]*leftShapeFunctions[leftIndex][col];
            for(CFuint i = 0; i < nbEqs; i++ )
            {
              elemMat(row*nbEqs + i, col*nbEqs + i)-=temp_value;
            }
           }
          else
          {

            CFreal temp_value=rightShapeFunctions[rightIndex][col-nbStatesInCellLeft]*rightShapeFunctions[rightIndex][row-nbStatesInCellLeft];
            for(CFuint i = 0; i < nbEqs; i++ )
            {
              elemMat(row*nbEqs + i, col*nbEqs + i)+=temp_value;
            }
           }
         }
      }
      elemMat*=leftWeight[0][kvadrature_point]*(m_diameter[leftID]/m_gOnMinusOne[leftID]*m_gRhoJump[leftID] + m_diameter[rightID]/m_gOnMinusOne[rightID]*m_gRhoJump[rightID])/2.0;//*detJacobi;
      acc.addValuesM(elemMat);

    }

    jacobMatrix->addValues(acc);

    // release the face
    geoBuilderFace->releaseGE();
  }


}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdStabilization::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_integrationIndex);
  result.push_back(&socket_normals);

  return result;
}

//////////////////////////////////////////////////////////////////////////////
  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

