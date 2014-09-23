#include "Framework/MethodCommandProvider.hh"
#include "Framework/CFSide.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"

#include "DiscontGalerkin/ViscousWallBC.hh"
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

MethodCommandProvider< ViscousWallBC,DiscontGalerkinSolverData,DiscontGalerkinModule >
  viscousWallBCProvider("ViscousWallBC");

//////////////////////////////////////////////////////////////////////////////

ViscousWallBC::ViscousWallBC(const std::string& name)
: ViscousBaseSolve(name),
    socket_rhs("rhs"),
    m_sigma(),
    m_theta()
{
}

//////////////////////////////////////////////////////////////////////////////

ViscousWallBC::~ViscousWallBC()
{
}

//////////////////////////////////////////////////////////////////////////////

void ViscousWallBC::setup()
{
  CFAUTOTRACE;
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  m_DF.resize(nbDim,nbEqs);
  // set a pointer to the cells
  m_cells.reset(MeshDataStack::getActive()->getTrs("InnerCells"));
  ViscousBaseSolve::setup();
  m_theta = getMethodData().getTheta();
  m_sigma = getMethodData().getSigma();
  // set Amatix to zero and set constant part of them
  m_kMatrix.resize(nbDim);
  for (CFuint k = 0; k < nbDim ; k++)
  {
    m_kMatrix[k].resize(nbDim);
    for (CFuint s = 0; s < nbDim ; s++)
    {
      m_kMatrix[k][s].resize(nbEqs,nbEqs);
      m_kMatrix[k][s]=0.;
    }
  }
  SafePtr<vector<ElementTypeData> > elementType = MeshDataStack::getActive()->getElementTypeData();

  // rhs storage
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

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

  m_mapElemData.sortKeys();

}

//////////////////////////////////////////////////////////////////////////////

void ViscousWallBC::compute_DFmatrix2D(State state, RealVector normal, RealMatrix *m_DF)
{
  CFreal v1=state[1]/state[0];
  CFreal v2=state[2]/state[0];

  (*m_DF)(0,0) = (v1*v1 + v2*v2)*normal[0]/2.0;
  (*m_DF)(0,1) = -v1*normal[0];
  (*m_DF)(0,2) = -v2*normal[0];
  (*m_DF)(0,3) = normal[0];

  (*m_DF)(1,0) = (v1*v1 + v2*v2)*normal[1]/2.0;
  (*m_DF)(1,1) = -v1*normal[1];
  (*m_DF)(1,2) = -v2*normal[1];
  (*m_DF)(1,3) = normal[1];

  return;
}

//////////////////////////////////////////////////////////////////////////////

void ViscousWallBC::compute_DFmatrix3D(State state, RealVector normal, RealMatrix *m_DF)
{
  CFreal v1=state[1]/state[0];
  CFreal v2=state[2]/state[0];
  CFreal v3=state[3]/state[0];

  (*m_DF)(0,0) = (v1*v1 + v2*v2 + v3*v3)*normal[0]/2.0;
  (*m_DF)(0,1) = -v1*normal[0];
  (*m_DF)(0,2) = -v2*normal[0];
  (*m_DF)(0,3) = -v3*normal[0];
  (*m_DF)(0,4) = normal[0];

  (*m_DF)(1,0) = (v1*v1 + v2*v2 + v3*v3)*normal[1]/2.0;
  (*m_DF)(1,1) = -v1*normal[1];
  (*m_DF)(1,2) = -v2*normal[1];
  (*m_DF)(1,3) = -v3*normal[1];
  (*m_DF)(1,4) = normal[1];

  (*m_DF)(2,0) = (v1*v1 + v2*v2 + v3*v3)*normal[2]/2.0;
  (*m_DF)(2,1) = -v1*normal[2];
  (*m_DF)(2,2) = -v2*normal[2];
  (*m_DF)(2,3) = -v3*normal[2];
  (*m_DF)(2,4) = normal[2];

  return;
}

//////////////////////////////////////////////////////////////////////////////

void ViscousWallBC::unsetup()
{
}
//////////////////////////////////////////////////////////////////////////////

void ViscousWallBC::executeOnTrs()
{
  CFAUTOTRACE;
  CFout << "ViscousWallBC applied to " << getCurrentTRS()->getName() << CFendl;
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  CFreal gamma=1.4;
  RealVector normal;
  State state;
  State Bstate;
  CFreal cv =1.0;

  SafePtr<LinearSystemSolver> lss = getMethodData().getLinearSystemSolver()[0];
  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();
  // rhs storage
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  Common::SafePtr<GeometricEntityPool<Framework::FaceToCellGEBuilder> >
     geoBuilder = getMethodData().getFaceBuilder();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get Wall Boundary Faces = CurrentTopologicalRegionSet
  SafePtr<TopologicalRegionSet> faces = getCurrentTRS();

  // get the geodata of the face builder and set the TRSs
  Framework::FaceToCellGEBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.cellsTRS = cells;
  geoData.facesTRS = faces;
  geoData.isBoundary = true;

  //get number of wall faces
  const CFuint nbFaces = faces->getLocalNbGeoEnts();
//   CFout << "  Faces:" << nbFaces << CFendl;

  //loop over all wall boundary faces
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace)
  {
//     CFout << "  Face:" << iFace << "   " << CFendl;
    //set index of face
    geoData.idx = iFace;
    //geo builder make face
    m_face = geoBuilder->buildGE();
    //get the (neighbouring) cell
    GeometricEntity* cellLeft  = m_face->getNeighborGeo(LEFT);
    //get nodes of face
    std::vector<Node*>&  nodes  = *m_face->getNodes();

    //get nodes and states of the cell
    std::vector<Node*>&  left_cell_nodes  = *cellLeft->getNodes();
    std::vector<State*>& left_cell_states = *cellLeft->getStates();

    if (left_cell_states[0]->isParUpdatable())
    {
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

      //get number of states of the cell
      const CFuint nbStatesInCellLeft = left_cell_states.size();

      //*****************************************************************
      //*****************************************************************
      //SET Integrator
      //compute shape function in quadrature points
      const std::vector<RealVector>& leftShapeFunctions =  getMethodData().getContourIntegrator()->getSolutionIntegrator(cellLeft)->computeShapeFunctionsAtQuadraturePoints();

      //numbers of quadrature points
      CFuint m_nbKvadrPoint = getMethodData().getContourIntegrator()->getSolutionIntegrator(cellLeft)->getIntegratorPattern()[0];

      //set weights for element quadrature
      const std::vector<RealVector>& leftWeight = getMethodData().getContourIntegrator()->getSolutionIntegrator(cellLeft)->getCoeff();

      //get coordinates of quadrature points
      const std::vector<RealVector>& leftCoord = getMethodData().getContourIntegrator()->getSolutionIntegrator(cellLeft)->getQuadraturePointsCoordinates();

      //compute gradient of shape functions in quadrature points
      std::vector<RealMatrix> leftGradient = cellLeft->computeSolutionShapeFunctionGradientsInMappedCoordinates(leftCoord);

      //*****************************************************************
      //*****************************************************************

      //get number of equations
      const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

      /// @todo this must be improved. Finding in a map
      /// for every cell is efficiently speaking not acceptable.
      //make block acumulator with size of states of the cell
      DGElemTypeData elemData = m_mapElemData.find(nbStatesInCellLeft);
      BlockAccumulator& acc = *elemData.first;
      RealMatrix& elemMat = *elemData.second;
      RealVector& elemVec = *elemData.third;

      //set matrix in blockaccumulator to 0
      elemMat=0.0;
      acc.setValuesM(elemMat);

      // set the IDs on the blockaccumulator (use setRowColIndex() )
      for (CFuint iState = 0; iState < nbStatesInCellLeft; ++iState) {
        const CFuint stateID = left_cell_states[iState]->getLocalID();
        acc.setRowColIndex(iState, stateID);
      }

      //find local index of face in cell (must be improved)
      CFuint m_idxFaceFromLeftCell= 10;
      if (nbDim == 2)
      {
        for (CFuint i=0; i < 3; i++)
          for (CFuint j=0; j < 2; j++)
            if((*left_cell_nodes[(i)%3]==*nodes[(j)%2])&&(*left_cell_nodes[(i+1)%3]==*nodes[(j+1)%2])) m_idxFaceFromLeftCell=i;
      }
      else
      {
        for (CFuint i=0; i < 4; i++)
          for (CFuint j=0; j < 3; j++)
            if(((*left_cell_nodes[(i)%4]==*nodes[(j)%3])&&(*left_cell_nodes[(i+1)%4]==*nodes[(j+1)%3])&&(*left_cell_nodes[(i+2)%4]==*nodes[(j+2)%3]))||((*left_cell_nodes[(i)%4]==*nodes[(j+2)%3])&&(*left_cell_nodes[(i+1)%4]==*nodes[(j+1)%3])&&(*left_cell_nodes[(i+2)%4]==*nodes[(j)%3])))
            {
              m_idxFaceFromLeftCell=i;
            }
      }
      cf_assert(m_idxFaceFromLeftCell!=10);

      //computation of normal (only for straight line connecting boundary nodes of face)
      normal.resize(nbDim);
      if (nbDim == 2)
      {
        normal[0]= (*left_cell_nodes[(m_idxFaceFromLeftCell+1) % 3])[YY] - (*left_cell_nodes[m_idxFaceFromLeftCell])[YY];
        normal[1]= (*left_cell_nodes[m_idxFaceFromLeftCell])[XX] - (*left_cell_nodes[(m_idxFaceFromLeftCell+1) % 3])[XX];
        normal /=sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
      }
      else
      {
        normal = m_face->computeAvgCellNormal();
        RealVector HlpNormal = *left_cell_nodes[(m_idxFaceFromLeftCell+3)%4] - *left_cell_nodes[(m_idxFaceFromLeftCell)%4];
        if ((normal[0]*HlpNormal[0]+normal[1]*HlpNormal[1]+normal[2]*HlpNormal[2]) > 0)
        {
          normal *=-1;
        }
        normal /=sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
      }
      //loop over kvadrature point on the face
      for(CFuint kvadrature_point = 0; kvadrature_point < m_nbKvadrPoint; kvadrature_point++ )
      {
        CFuint leftIndex  = m_idxFaceFromLeftCell*m_nbKvadrPoint + kvadrature_point;
        elemMat=0.0;
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) //loop over members of state - set to zero
        {
          state[iEq]= 0.0;
          Bstate[iEq]= 0.0;
        }
        //computation of state in point of kvadrature - from previous time step
        for (CFuint iState = 0; iState < nbStatesInCellLeft; ++iState) //loop over states in cell
        {
          RealVector &states = *left_cell_states[iState]->getData();
          for (CFuint iEq = 0; iEq < nbEqs; ++iEq) //loop over members of state
          {
            state[iEq] += leftShapeFunctions[leftIndex][iState]*states[iEq];
          }
        }
        if (nbDim == 2)
        {
          compute_Kmatrix2D(state,&(m_kMatrix));
        }
        else
        {
          compute_Kmatrix3D(state,&(m_kMatrix));
        }
        Bstate[0] = state[0];
        for (CFuint iEq = 1; iEq < nbEqs-1; ++iEq) //loop over members of state
          {
            Bstate[nbEqs-1] -= state[iEq]*state[iEq]/(state[0]);
          }
        Bstate[nbEqs-1] = (Bstate[nbEqs-1]/2.0 + state[nbEqs-1])/cv;
        //compute DF matrix from numerical flux - wall term
        if (nbDim == 2)
        {
          compute_DFmatrix2D(state, normal, &m_DF);
        }
        else
        {
          compute_DFmatrix3D(state, normal, &m_DF);
        }
        //loop over test function
        for(CFuint row = 0; row < nbStatesInCellLeft; row++ )
        {
          //loop over base function of solution
          for(CFuint col = 0; col < nbStatesInCellLeft; col++ )
          {
            //temp_value = multiplication of test function and base function in kvadrature point
            CFreal temp_value=leftShapeFunctions[leftIndex][row]*leftShapeFunctions[leftIndex][col];
            for(CFuint i = 0; i < nbEqs-2; i++ )
              for(CFuint j = 0; j < nbEqs; j++ )
              {
                elemMat(row*nbEqs + i + 1, col*nbEqs + j)+=(m_DF(i,j))*temp_value;
              }
          }
        } // end of numerical flux
        elemMat*=leftWeight[0][kvadrature_point]*(gamma - 1)*detJacobi;
        acc.addValuesM(elemMat);

        elemMat=0.0;
        if (kvadrature_point == 0) elemVec=0.0;
        for(CFuint row = 0; row < nbStatesInCellLeft; row++ )
        {
//         loop over base function of solution
          for(CFuint col = 0; col < nbStatesInCellLeft; col++ )
          {
            for(CFuint k=0;k<nbDim;k++)
            {
              for(CFuint s=0;s<nbDim;s++)
              {
//             temp_value = multiplication of test function, base function in kvadrature point and normal
              CFreal temp_value=leftGradient[leftIndex](col,k)*leftShapeFunctions[leftIndex][row]*normal[s];
              for(CFuint i = 0; i < nbEqs; i++ )
                for(CFuint j = 0; j < nbEqs; j++ )
                {
//                 K_{sk} * \frac{(\partial w}{\partial x_k}|_L * \varphi|_L
                  elemMat(row*nbEqs + i, col*nbEqs + j)+=m_kMatrix[s][k](i,j)*temp_value;
//                  m_theta * K_{sk} * \frac{(\partial \varphi}{\partial x_k}|_L * \w|_L
                  elemMat(col*nbEqs + j, row*nbEqs + i)+=m_theta*m_kMatrix[s][k](i,j)*temp_value;
                }
              }
            }
          }
          //add part from local rhs to global rhs
          const CFuint stateID = left_cell_states[row]->getLocalID();
          for(CFuint k=0;k<nbDim;k++)
          {
            for(CFuint s=0;s<nbDim;s++)
            {
              for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
              {
                for (CFuint j = 0; j < nbEqs; ++j)
                {
                  rhs(stateID, iEq, nbEqs) -= m_theta*m_kMatrix[s][k](j,iEq)*leftGradient[leftIndex](row,k)*Bstate[j]*normal[s]*leftWeight[0][kvadrature_point]*detJacobi;
                }
              }
            }
          }
          for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
          {
            rhs(stateID, iEq, nbEqs) += m_sigma*leftShapeFunctions[leftIndex][row]*Bstate[iEq]*leftWeight[0][kvadrature_point]*detJacobi/detJacobi/m_Re; //Cw/gamma
//             elemVec[row + iEq*nbStatesInCellLeft] += m_sigma*leftShapeFunctions[leftIndex][row]*Bstate[iEq]*leftWeight[0][kvadrature_point]*detJacobi/detJacobi/m_Re; //Cw/gamma
          }
          for(CFuint col = 0; col < nbStatesInCellLeft; col++ )
          {
            for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
            {
              elemMat(row*nbEqs + iEq, col*nbEqs + iEq)-=m_sigma*leftShapeFunctions[leftIndex][row]*leftShapeFunctions[leftIndex][col]/detJacobi/m_Re; //Cw/gamma
            }
          }
        }
        elemMat*=-leftWeight[0][kvadrature_point]*detJacobi;
        acc.addValuesM(elemMat);
//         acc.printToScreen();
      }
      // add the values in the jacobian matrix
//       if (iFace == 6)
//          CFout << "\n" << elemVec << "\n";
//         acc.printToScreen();
      jacobMatrix->addValues(acc);
    }
    // release the face
    geoBuilder->releaseGE();
  }
//   jacobMatrix->finalAssembly();
//   jacobMatrix->printToFile("inside");

  CFout << " ... OK \n" << CFendl;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> > ViscousWallBC::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

