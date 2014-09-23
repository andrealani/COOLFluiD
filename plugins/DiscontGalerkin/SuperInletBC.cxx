#include "Framework/MethodCommandProvider.hh"
#include "Framework/CFSide.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/FaceToCellGEBuilder.hh"

#include "DiscontGalerkin/SuperInletBC.hh"
#include "DiscontGalerkin/DiscontGalerkin.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< SuperInletBC,DiscontGalerkinSolverData,DiscontGalerkinModule >
  SuperInletBCProvider("SuperInletBC");

//////////////////////////////////////////////////////////////////////////////

void SuperInletBC::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
}

//////////////////////////////////////////////////////////////////////////////

SuperInletBC::SuperInletBC(const std::string& name)
: StdBaseSolve(name),
  socket_rhs("rhs"),
  kappa1(0.4)

{
  addConfigOptionsTo(this);
  _functions = std::vector<std::string>();
  setParameter("Def",&_functions);

  _vars = std::vector<std::string>();
  setParameter("Vars",&_vars);
}

//////////////////////////////////////////////////////////////////////////////

SuperInletBC::~SuperInletBC()
{
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletBC::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DiscontGalerkinSolverCom::configure(args);

  _vFunction.setFunctions(_functions);
  _vFunction.setVariables(_vars);
  try {
    _vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletBC::setup()
{
  CFAUTOTRACE;

  // set a pointer to the cells
  m_cells.reset(MeshDataStack::getActive()->getTrs("InnerCells"));

  SafePtr<vector<ElementTypeData> > elementType = MeshDataStack::getActive()->getElementTypeData();

  // rhs storage
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  const CFuint nbElemTypes = elementType->size();
   const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

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


void SuperInletBC::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void SuperInletBC::executeOnTrs()
{
  CFAUTOTRACE;
  CFout << "SuperInletBC applied to " << getCurrentTRS()->getName() << CFendl;
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  RealMatrix T(nbEqs,nbEqs);
  RealMatrix T1(nbEqs,nbEqs);
  RealMatrix Pplus(nbEqs,nbEqs);
  RealMatrix Pminus(nbEqs,nbEqs);
  RealMatrix EigenVal(2,nbEqs);
  RealVector normal;
  State  state;
  State  bstate;
  State  Dstate;
  Node   Dnode;

  SafePtr<LinearSystemSolver> lss = getMethodData().getLinearSystemSolver()[0];
  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  Common::SafePtr<GeometricEntityPool<Framework::FaceToCellGEBuilder> >
     geoBuilder = getMethodData().getFaceBuilder();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get Inlet Boundary Faces = CurrentTopologicalRegionSet
  SafePtr<TopologicalRegionSet> faces = getCurrentTRS();

  // get the geodata of the face builder and set the TRSs
  Framework::FaceToCellGEBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.cellsTRS = cells;
  geoData.facesTRS = faces;
  geoData.isBoundary = true;

  //get number of inlet faces
  const CFuint nbFaces = faces->getLocalNbGeoEnts();
//   CFout << "  Faces:" << nbFaces << CFendl;

  //loop over all inlet boundary faces
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace)
  {
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

    //computation of the inverse Jacobi matrix of mapping from refference element to the cell
    //invJacobiLeft[0] = (*left_cell_nodes[2])[YY] - (*left_cell_nodes[0])[YY];
    //invJacobiLeft[1] = - (*left_cell_nodes[2])[XX] + (*left_cell_nodes[0])[XX];
    //invJacobiLeft[2] = - (*left_cell_nodes[1])[YY] + (*left_cell_nodes[0])[YY];
    //invJacobiLeft[3] = (*left_cell_nodes[1])[XX] - (*left_cell_nodes[0])[XX];

    //computation of the Jacobi determinant of mapping from refference element to cell
    //detJacobiLeft = invJacobiLeft[0]*invJacobiLeft[3] - invJacobiLeft[1]*invJacobiLeft[2];
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
    }
    else
    {

      normal = m_face->computeAvgCellNormal();

      RealVector HlpNormal = *left_cell_nodes[(m_idxFaceFromLeftCell+3)%4] - *left_cell_nodes[(m_idxFaceFromLeftCell)%4];

      if ((normal[0]*HlpNormal[0]+normal[1]*HlpNormal[1]+normal[2]*HlpNormal[2]) > 0)
      {
        normal *=-1;
      }
    }

    //loop over kvadrature point on the face
    for(CFuint kvadrature_point = 0; kvadrature_point < m_nbKvadrPoint; kvadrature_point++ )
    {
      CFuint leftIndex  = m_idxFaceFromLeftCell*m_nbKvadrPoint + kvadrature_point;
      elemMat=0.0;
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) //loop over members of state - set to zero
      {
        state[iEq]= 0;
      }
      //computation of state in point of kvadrature - from previous time step

      for (CFuint iState = 0; iState < nbStatesInCellLeft; ++iState) //loop over states in cell
      {
        RealVector &states = *left_cell_states[iState]->getData();
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) //loop over members of state
        {
          state[iEq] += leftShapeFunctions[kvadrature_point+m_idxFaceFromLeftCell*m_nbKvadrPoint][iState]*states[iEq];
        }
      }
      //set coordinates of nodes for gauss kvadrature - to evaluate function from boundary condition

      if (nbDim == 2)
      {
        Dnode = leftCoord[m_idxFaceFromLeftCell*m_nbKvadrPoint+kvadrature_point][XX]*(*left_cell_nodes[0])
             +leftCoord[m_idxFaceFromLeftCell*m_nbKvadrPoint+kvadrature_point][YY]*(*left_cell_nodes[1])
             +(1-leftCoord[m_idxFaceFromLeftCell*m_nbKvadrPoint+kvadrature_point][XX]-leftCoord[m_idxFaceFromLeftCell*m_nbKvadrPoint+kvadrature_point][YY])*(*left_cell_nodes[2]);
      }
      else
      {
        Dnode = leftCoord[m_idxFaceFromLeftCell*m_nbKvadrPoint+kvadrature_point][XX]*(*left_cell_nodes[0])
             +leftCoord[m_idxFaceFromLeftCell*m_nbKvadrPoint+kvadrature_point][YY]*(*left_cell_nodes[1])
             +leftCoord[m_idxFaceFromLeftCell*m_nbKvadrPoint+kvadrature_point][ZZ]*(*left_cell_nodes[2])
             +(1-leftCoord[m_idxFaceFromLeftCell*m_nbKvadrPoint+kvadrature_point][XX]   -leftCoord[m_idxFaceFromLeftCell*m_nbKvadrPoint+kvadrature_point][YY] -leftCoord[m_idxFaceFromLeftCell*m_nbKvadrPoint+kvadrature_point][ZZ])*(*left_cell_nodes[3]);
      }

      //Dstate = evaluate function from boundary condition in Dnode coordinates
      _vFunction.evaluate(Dnode,Dstate);

      //compute matrixes P+ a P- in point of kvadrature
      if (nbDim == 2)
      {
        compute_EigenValVec2D(state, T, T1, &Pplus, &Pminus, &EigenVal, normal);
      }
      else
      {
        if (compute_EigenValVec3D(state, T, T1, &Pplus, &Pminus, &EigenVal, normal) !=0)
        {
          Node Dnode;
          Dnode = leftCoord[leftIndex][0]*(*(left_cell_nodes[0]))
               +leftCoord[leftIndex][1]*(*(left_cell_nodes[1]))
               +leftCoord[leftIndex][2]*(*(left_cell_nodes[2]))
               +(1-leftCoord[leftIndex][0]-leftCoord[leftIndex][1] -leftCoord[leftIndex][2])*(*(left_cell_nodes[3]));
          CFout << "\n  Negative pressure in " << Dnode << "  STATE  "  << state << "  normal  " << normal << CFendl;
        }
      }
      //compute inner boundary face term 
      //loop over test function
      for(CFuint row = 0; row < nbStatesInCellLeft; row++ )
      {
        //loop over base function of solution
        for(CFuint col = 0; col < nbStatesInCellLeft; col++ )
        {
          //temp_value = multiplication of test function and base function in kvadrature point
          CFreal temp_value=leftShapeFunctions[kvadrature_point+m_idxFaceFromLeftCell*m_nbKvadrPoint][row]*leftShapeFunctions[kvadrature_point+m_idxFaceFromLeftCell*m_nbKvadrPoint][col];
          for(CFuint i = 0; i < nbEqs; i++ )
            for(CFuint j = 0; j < nbEqs; j++ )
            {
              elemMat(row*nbEqs + i, col*nbEqs + j)+=(Pplus(i,j))*temp_value;
            }
        }
      }
      elemMat*=leftWeight[0][kvadrature_point]*detJacobi;
      acc.addValuesM(elemMat);

      elemVec = 0.0;
      //compute inner boundary face term -> to RHS
      //loop over test function
      for(CFuint row = 0; row < nbStatesInCellLeft; row++ )
      {
        for(CFuint i = 0; i < nbEqs; i++ )//loop over member of state - to test function
        {
          elemVec[row*nbEqs+i]=0; //set RHS[index]
          for(CFuint j = 0; j < nbEqs; j++ )//loop over member of stateBC
          {
            elemVec[row*nbEqs+i]+=(Pminus(i,j))*(Dstate)[j];
          }
          elemVec[row*nbEqs+i]*=leftShapeFunctions[kvadrature_point+m_idxFaceFromLeftCell*m_nbKvadrPoint][row];
        }
      } // end of numerical flux
      for (CFuint iState = 0; iState < nbStatesInCellLeft; ++iState) {
        const CFuint stateID = left_cell_states[iState]->getLocalID();
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
          rhs(stateID, iEq, nbEqs) -= elemVec[iState*nbEqs + iEq]*leftWeight[0][kvadrature_point]*detJacobi;
// CFout << "  EL:  " << rhs(stateID, iEq, nbEqs) << CFendl;
        }
      }
    }

    // add the values in the jacobian matrix
    jacobMatrix->addValues(acc);
    // release the face
    geoBuilder->releaseGE();
  }
//   jacobMatrix->finalAssembly();

  CFout << " ... OK\n" << CFendl;
//       jacobMatrix->finalAssembly();   
//       jacobMatrix->printToFile("inside");
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> > SuperInletBC::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

