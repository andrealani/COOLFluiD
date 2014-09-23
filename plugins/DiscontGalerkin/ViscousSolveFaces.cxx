#include "Framework/MethodCommandProvider.hh"
#include "Framework/CFSide.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"

#include "DiscontGalerkin/ViscousSolveFaces.hh"
#include "DiscontGalerkin/DiscontGalerkin.hh"
#include "Framework/FaceToCellGEBuilder.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

namespace COOLFluiD {
  namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ViscousSolveFaces,DiscontGalerkinSolverData,DiscontGalerkinModule >
  viscousSolveFacesProvider("ViscousSolveFaces");

//////////////////////////////////////////////////////////////////////////////

ViscousSolveFaces::ViscousSolveFaces(const std::string& name)
: ViscousBaseSolve(name),
    socket_rhs("rhs"),
    socket_integrationIndex("integrationIndex"),
    socket_normals("normals")
{
}

//////////////////////////////////////////////////////////////////////////////

ViscousSolveFaces::~ViscousSolveFaces()
{
}

//////////////////////////////////////////////////////////////////////////////

void ViscousSolveFaces::setup()
{
  CFAUTOTRACE;
  ViscousBaseSolve::setup();
  m_Theta = getMethodData().getTheta();
  m_Sigma = getMethodData().getSigma();
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  // set Kmatix to zero
  m_kMatrix.resize(2);
  for (CFuint i = 0; i < 2 ; i++)
  {
    m_kMatrix[i].resize(nbDim);
    for (CFuint k = 0; k < nbDim ; k++)
    {
      m_kMatrix[i][k].resize(nbDim);
      for (CFuint s = 0; s < nbDim ; s++)
      {
        m_kMatrix[i][k][s].resize(nbEqs,nbEqs);
        m_kMatrix[i][k][s]=0.;
      }
    }
  }
//   CFreal gamma = 1.4;

  // set a pointer to the cells
  m_cells.reset(MeshDataStack::getActive()->getTrs("InnerCells"));

  //get types of elements in triangulation
  SafePtr<vector<ElementTypeData> > elementType = MeshDataStack::getActive()->getElementTypeData();

  // rhs storage
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  const CFuint nbElemTypes = elementType->size();

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

}

//////////////////////////////////////////////////////////////////////////////

void ViscousSolveFaces::unsetup()
{
  CFAUTOTRACE;
  ViscousBaseSolve::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void ViscousSolveFaces::execute()
{
  CFAUTOTRACE;
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  RealMatrix T(nbEqs,nbEqs);
  RealMatrix T1(nbEqs,nbEqs);
  RealMatrix Pplus(nbEqs,nbEqs);
  RealMatrix Pminus(nbEqs,nbEqs);
  RealMatrix EigenVal(2,nbEqs);
  RealVector normal;
  State  stateA;

  DataHandle< std::vector< CFuint > >
    integrationIndex = socket_integrationIndex.getDataHandle();
  DataHandle< RealVector > normals = socket_normals.getDataHandle();
  SafePtr<LinearSystemSolver> lss = getMethodData().getLinearSystemSolver()[0];
  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();
  // reset to zero all non zero entries in the jacobian matrix
  jacobMatrix->resetToZeroEntries();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  Common::SafePtr<GeometricEntityPool<Framework::FaceToCellGEBuilder> >
     geoBuilder = getMethodData().getFaceBuilder();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get InnerFaces TopologicalRegionSet
  SafePtr<TopologicalRegionSet> faces = MeshDataStack::getActive()->getTrs("InnerFaces");

  CFout << "ViscousSolveFaces applied to " << faces->getName() << CFendl;

  // get the geodata of the face builder and set the TRSs
  Framework::FaceToCellGEBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.cellsTRS = cells;
  geoData.facesTRS = faces;
  geoData.isBoundary = false;

  const CFuint nbFaces = faces->getLocalNbGeoEnts();
  //loop over all inner faces

  for (CFuint iFace = 0; iFace < nbFaces; ++iFace)
  {
    CFLogDebugMax("Face " << iFace << "\n");
    //set index of face
    geoData.idx = iFace;
    //geo builder make face
    GeometricEntity& face = *geoBuilder->buildGE();

    //nodes of the actual face
    std::vector<Node*>&  nodes  = *face.getNodes();
    //get the left (neighbouring) cell
    GeometricEntity* cellLeft  = face.getNeighborGeo(LEFT);

//     if ((nodes[0]->getLocalID() == 16) || (nodes[0]->getLocalID() == 28))
//     CFout << "\n" << iFace << "  " << nodes[0]->getLocalID() << "  " << nodes[1]->getLocalID();
    //get nodes and states of the left cell
    std::vector<Node*>&  left_cell_nodes  = *cellLeft->getNodes();
    std::vector<State*>& left_cell_states = *cellLeft->getStates();

    //get the right (neighbouring) cell
    GeometricEntity* cellRight = face.getNeighborGeo(RIGHT);
    //get nodes and states of the right cell
//     std::vector<Node*>&  right_cell_nodes  = *cellRight->getNodes();
    std::vector<State*>& right_cell_states = *cellRight->getStates();

//     get number of states in left and right cell
    const CFuint nbStatesInCellLeft  = left_cell_states.size();
    const CFuint nbStatesInCellRight = right_cell_states.size();

    CFreal avg_massCell = (abs(cellRight->computeVolume()) + abs(cellLeft->computeVolume()))/2.0;
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
//     RealVector& elemVec = *elemData.third;

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
      cf_assert((hlpVector.norm1()/detJacobi < 0.0001));

      //LEFT CELL
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) //loop over members of state - set to zero
      {
        (*m_state)[iEq] = 0;
      }
      //computation of state in point of kvadrature - from previous time step
      for (CFuint iState = 0; iState < nbStatesInCellLeft; ++iState) //loop over states in cell
      {
        RealVector &states = *left_cell_states[iState]->getData();
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) //loop over members of state
        {
          (*m_state)[iEq] += leftShapeFunctions[leftIndex][iState]*states[iEq];
        }
      }
      //call compute Kmatrix from baseSolve for LEFT cell
      if (nbDim == 2)
      {
        compute_Kmatrix2D(*m_state,&(m_kMatrix[LEFT]));
      }
      else
      {
        compute_Kmatrix3D(*m_state,&(m_kMatrix[LEFT]));
      }
      //half of computation of average of state in point of kvadrature
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) //loop over members of state - set to zero
      {
        (stateA)[iEq] = (*m_state)[iEq];
      }

      //RIGHT CELL
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) //loop over members of state
      {
        (*m_state)[iEq] = 0;
      }
      //computation of state and gradient of state in point of kvadrature - from previous step
      for (CFuint iState = 0; iState < nbStatesInCellRight; ++iState) //loop over states in cell 
      {
        RealVector &states = *right_cell_states[iState]->getData();
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) //loop over members of state
        {
          (*m_state)[iEq] += rightShapeFunctions[rightIndex][iState]*states[iEq];
        }
      }
      //call compute Kmatrix from baseSolve for RIGHT cell
      if (nbDim == 2)
      {
        compute_Kmatrix2D(*m_state,&(m_kMatrix[RIGHT]));
      }
      else
      {
        compute_Kmatrix3D(*m_state,&(m_kMatrix[RIGHT]));
      }
      //second half of computation of average of state in point of kvadrature
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) //loop over members of state
      {
        (stateA)[iEq] += (*m_state)[iEq];
        (stateA)[iEq]/=2.0;
      }

      // part for viscous term
      CFreal temp_value;
      elemMat=0.0;
      // LL
      for(CFuint col = 0; col < nbStatesInCellLeft; col++ )
      {
//         loop over base function of solution
        for(CFuint row = 0; row < nbStatesInCellLeft; row++ )
        {
          for(CFuint k=0;k<nbDim;k++)
          {
            for(CFuint s=0;s<nbDim;s++)
            {
//             temp_value = multiplication of test function, base function in kvadrature point and normal
            temp_value=leftGradient[leftIndex](col,k)*leftShapeFunctions[leftIndex][row]*normal[s];
            for(CFuint i = 0; i < nbEqs; i++ )
              for(CFuint j = 0; j < nbEqs; j++ )
              {
//                 K_{sk} * \frac{(\partial w}{\partial x_k}|_L * \varphi|_L
                elemMat(row*nbEqs + i, col*nbEqs + j)+=m_kMatrix[LEFT][s][k](i,j)*temp_value;
//                  m_Theta * K_{sk} * \frac{(\partial \varphi}{\partial x_k}|_L * \w|_L
                elemMat(col*nbEqs + j, row*nbEqs + i)+=m_Theta*m_kMatrix[LEFT][s][k](i,j)*temp_value;
              }
            }
          }
        }
      // LR
        for(CFuint row = nbStatesInCellLeft; row < nbStatesInCellLeft+nbStatesInCellRight; row++ )
        {
          for(CFuint k=0;k<nbDim;k++)
          {
            for(CFuint s=0;s<nbDim;s++)
            {
//             temp_value = multiplication of test function, base function in kvadrature point and normal
            temp_value=leftGradient[leftIndex](col,k)*rightShapeFunctions[rightIndex][row-nbStatesInCellLeft]*normal[s];
            for(CFuint i = 0; i < nbEqs; i++ )
              for(CFuint j = 0; j < nbEqs; j++ )
              {
//                  K_{sk} * \frac{(\partial w}{\partial x_k}|_L * \varphi|_P
                elemMat(row*nbEqs + i, col*nbEqs + j)-=m_kMatrix[LEFT][s][k](i,j)*temp_value;
//                  m_Theta * K_{sk} * \frac{(\partial \varphi}{\partial x_k}|_L * \w|_P
                elemMat(col*nbEqs + j, row*nbEqs + i)-=m_Theta*m_kMatrix[RIGHT][s][k](i,j)*temp_value;
              }
            }
          }
        }
      }
      // RL
      for(CFuint col = nbStatesInCellLeft; col < nbStatesInCellLeft + nbStatesInCellRight; col++ )
      {
//         loop over base function of solution
        for(CFuint row = 0; row < nbStatesInCellLeft; row++ )
        {
          for(CFuint k=0;k<nbDim;k++)
          {
            for(CFuint s=0;s<nbDim;s++)
            {
//             temp_value = multiplication of test function, base function in kvadrature point and normal
            temp_value=rightGradient[rightIndex](col-nbStatesInCellLeft,k)*leftShapeFunctions[leftIndex][row]*normal[s];
            for(CFuint i = 0; i < nbEqs; i++ )
              for(CFuint j = 0; j < nbEqs; j++ )
              {
//                  K_{sk} * \frac{(\partial w}{\partial x_k}|_P * \varphi|_L
                elemMat(row*nbEqs + i, col*nbEqs + j)+=m_kMatrix[RIGHT][s][k](i,j)*temp_value;
//                  m_Theta * K_{sk} * \frac{(\partial \varphi}{\partial x_k}|_P * \w|_L
                elemMat(col*nbEqs + j, row*nbEqs + i)+=m_Theta*m_kMatrix[LEFT][s][k](i,j)*temp_value;
              }
            }
          }
        }
        //RR
        for(CFuint row = nbStatesInCellLeft; row < nbStatesInCellLeft+nbStatesInCellRight; row++ )
        {
          for(CFuint k=0;k<nbDim;k++)
          {
            for(CFuint s=0;s<nbDim;s++)
            {
//             temp_value = multiplication of test function, base function in kvadrature point and normal
            temp_value=rightGradient[rightIndex](col-nbStatesInCellLeft,k)*rightShapeFunctions[rightIndex][row-nbStatesInCellLeft]*normal[s];
            for(CFuint i = 0; i < nbEqs; i++ )
              for(CFuint j = 0; j < nbEqs; j++ )
              {
//                  K_{sk} * \frac{(\partial w}{\partial x_k}|_P * \varphi|_P
                elemMat(row*nbEqs + i, col*nbEqs + j)-=m_kMatrix[RIGHT][s][k](i,j)*temp_value;
//                  m_Theta * K_{sk} * \frac{(\partial \varphi}{\partial x_k}|_P * \w|_P
                elemMat(col*nbEqs + j, row*nbEqs + i)-=m_Theta*m_kMatrix[RIGHT][s][k](i,j)*temp_value;
              }
            }
          }
        }
      }
      elemMat*=-0.5*leftWeight[0][kvadrature_point]*detJacobi;
      acc.addValuesM(elemMat);

      //set element matrix to zero
      elemMat=0.0;
      //compute matrixes P+ a P- in point of kvadrature
      if (nbDim == 2)
      {
        compute_EigenValVec2D(stateA, T, T1, &Pplus, &Pminus, &EigenVal, normal);
        if (getMethodData().getMaxEigenval() < abs(EigenVal(0,2))*detJacobi/avg_massCell)
        {
          getMethodData().setMaxEigenval(abs(EigenVal(0,2))*detJacobi/avg_massCell);
        }
      }
      else
      {
        if (compute_EigenValVec3D(stateA, T, T1, &Pplus, &Pminus, &EigenVal, normal) !=0)
        {
          Node Dnode;
          Dnode = leftCoord[leftIndex][0]*(*(left_cell_nodes[0]))
               +leftCoord[leftIndex][1]*(*(left_cell_nodes[1]))
               +leftCoord[leftIndex][2]*(*(left_cell_nodes[2]))
               +(1-leftCoord[leftIndex][0]-leftCoord[leftIndex][1] -leftCoord[leftIndex][2])*(*(left_cell_nodes[3]));
          CFout << "\n  Negative pressure in " << Dnode << "  STATE  "  << stateA << "  normal  " << normal << CFendl;
        }
        if (getMethodData().getMaxEigenval() < abs(EigenVal(0,4))*detJacobi/avg_massCell)
        {
          getMethodData().setMaxEigenval(abs(EigenVal(0,4))*detJacobi/avg_massCell);
        }
      }
      //compute inner face term
      //loop over test function
//       CFreal temp_value;
      for(CFuint row = 0; row < nbStatesInCellLeft; row++ )
      {
        //loop over base function of solution
        for(CFuint col = 0; col < nbStatesInCellLeft + nbStatesInCellRight; col++ )
        {
          //test and base functions are from left cell
          if (col < nbStatesInCellLeft)
          {
if (left_cell_states[0]->isParUpdatable())
{

            //temp_value = multiplication of test function and base function in kvadrature point
            temp_value=leftShapeFunctions[leftIndex][col]*leftShapeFunctions[leftIndex][row];
            for(CFuint i = 0; i < nbEqs; i++ )
            {
              for(CFuint j = 0; j < nbEqs; j++ )
              {
                elemMat(row*nbEqs + i, col*nbEqs + j)+=(Pplus(i,j))*temp_value;
              }
              elemMat(row*nbEqs + i, col*nbEqs + i)+=m_Sigma*temp_value/detJacobi/m_Re; //Cw/gamma 
            }
}
          }
          //test function is from left cell and base functin from right cell
          else
          {
if (right_cell_states[0]->isParUpdatable())
{
            temp_value=rightShapeFunctions[rightIndex][col-nbStatesInCellLeft]*leftShapeFunctions[leftIndex][row];
            for(CFuint i = 0; i < nbEqs; i++ )
            {
              for(CFuint j = 0; j < nbEqs; j++ )
              {
                elemMat(row*nbEqs + i, col*nbEqs + j)+=(Pminus(i,j))*temp_value;
              }
              elemMat(row*nbEqs + i, col*nbEqs + i)-=m_Sigma*temp_value/detJacobi/m_Re; //Cw/gamma 
            }
}
          }
        }
      }
//       normal *=-1;
//       if (nbDim == 2)
//       {
//         compute_EigenValVec2D(stateA, T, T1, &Pplus, &Pminus, &EigenVal, normal);
//         if (getMethodData().getMaxEigenval() < abs(EigenVal(0,2))*detJacobi/avg_massCell)
//         {
//           getMethodData().setMaxEigenval(abs(EigenVal(0,2))*detJacobi/avg_massCell);
//         }
//       }
//       else
//       {
//          compute_EigenValVec3D(stateA, T, T1, &Pplus, &Pminus, &EigenVal, normal);
//         if (getMethodData().getMaxEigenval() < abs(EigenVal(0,4))*detJacobi/avg_massCell)
//         {
//           getMethodData().setMaxEigenval(abs(EigenVal(0,4))*detJacobi/avg_massCell);
//         }
//       }
      for(CFuint row = nbStatesInCellLeft; row < nbStatesInCellLeft + nbStatesInCellRight; row++ )
      {
        //loop over base function of solution
        for(CFuint col = 0; col < nbStatesInCellLeft + nbStatesInCellRight; col++ )
        {
          //test function is from right cell and base functin from left cell
          if (col < nbStatesInCellLeft)
          {
if (left_cell_states[0]->isParUpdatable())
{
            temp_value=leftShapeFunctions[leftIndex][col]*rightShapeFunctions[rightIndex][row-nbStatesInCellLeft];
            for(CFuint i = 0; i < nbEqs; i++ )
            {
              for(CFuint j = 0; j < nbEqs; j++ )
              {
//                 elemMat(row*nbEqs + i, col*nbEqs + j)+=(Pminus(i,j))*temp_value;
                elemMat(row*nbEqs + i, col*nbEqs + j)-=(Pplus(i,j))*temp_value;
              }
              elemMat(row*nbEqs + i, col*nbEqs + i)-=m_Sigma*temp_value/detJacobi/m_Re; //Cw/gamma 
            }
}
          }
          //test and base functions are from right cell
          else
          {
if (right_cell_states[0]->isParUpdatable())
{
            temp_value=rightShapeFunctions[rightIndex][col-nbStatesInCellLeft]*rightShapeFunctions[rightIndex][row-nbStatesInCellLeft];
            for(CFuint i = 0; i < nbEqs; i++ )
            {
              for(CFuint j = 0; j < nbEqs; j++ )
              {
//                 elemMat(row*nbEqs + i, col*nbEqs + j)+=(Pplus(i,j))*temp_value;
                elemMat(row*nbEqs + i, col*nbEqs + j)-=(Pminus(i,j))*temp_value;
              }
              elemMat(row*nbEqs + i, col*nbEqs + i)+=m_Sigma*temp_value/detJacobi/m_Re; //Cw/gamma 
            }
}
          }
        }
      } // end of numerical flux
//       normal *=-1;
      elemMat*=leftWeight[0][kvadrature_point]*detJacobi;
      acc.addValuesM(elemMat);
    }
// if (iFace == 185)
//    acc.printToScreen();
    // add the values in the jacobian matrix
    jacobMatrix->addValues(acc);
    // release the face
    geoBuilder->releaseGE();
  }
//   jacobMatrix->finalAssembly();
//   jacobMatrix->printToFile("inside");
CFout <<  " ... OK\n" << CFendl;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ViscousSolveFaces::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_integrationIndex);
  result.push_back(&socket_normals);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

