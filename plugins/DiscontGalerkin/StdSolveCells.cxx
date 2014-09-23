#include "Framework/MethodCommandProvider.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/SubSystemStatus.hh"

#include "DiscontGalerkin/StdSolveCells.hh"
#include "DiscontGalerkin/DiscontGalerkin.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< StdSolveCells,DiscontGalerkinSolverData,DiscontGalerkinModule >
  stdSolveCellsProvider("StdSolveCells");

//////////////////////////////////////////////////////////////////////////////

StdSolveCells::StdSolveCells(const std::string& name)
  : StdBaseSolve(name),
    m_mapElemData(),
    socket_rhs("rhs"),
    socket_states("states"),
    socket_old_states("old_states"),
    m_cells()
{
}

//////////////////////////////////////////////////////////////////////////////

StdSolveCells::~StdSolveCells()
{
}

//////////////////////////////////////////////////////////////////////////////

void StdSolveCells::setup()
{
  CFAUTOTRACE;
  StdBaseSolve::setup();
  m_Alpha = getMethodData().getAlpha();
  m_MaxCFL = getMethodData().getMaxCFL();
//   CFreal gamma = 1.4;
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  // set Amatix to zero and set constant part of them
  m_aMatrix.resize(nbDim);
  for (CFuint i =0; i<nbDim; i++)
  {
    m_aMatrix[i].resize(nbEqs,nbEqs);
    m_aMatrix[i] = 0.;
  }

  // set a pointer to the inner cells
  m_cells.reset(MeshDataStack::getActive()->getTrs("InnerCells"));

  //get types of elements in triangulation
  SafePtr<vector<ElementTypeData> > elementType = MeshDataStack::getActive()->getElementTypeData();

  // rhs storage
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  m_oldState = new Framework::State();

  DataHandle<State*> old_states = socket_old_states.getDataHandle();
  const CFuint nbstates = rhs.size()/nbEqs;
  old_states.resize(nbstates);
  for(CFuint i=0; i<nbstates; ++i)
  {
    old_states[i] = new State();
  }

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

void StdSolveCells::unsetup()
{
  CFAUTOTRACE;
  StdBaseSolve::unsetup();

  // deallocate our memory
  DataHandle<State*> old_states = socket_old_states.getDataHandle();
  const CFuint nbstates = old_states.size();
  for(CFuint i=0; i < nbstates; ++i)
  {
    deletePtr ( old_states[i] );
  }
}
//////////////////////////////////////////////////////////////////////////////

void StdSolveCells::execute()
{
  CFAUTOTRACE;
  CFuint m_nbKvadrPoint;
  CFout << "StdSolveCells applied to " << m_cells->getName() << CFendl;
  static CFreal tau = 0.0;
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();

  // get rhs
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  DataHandle<State*> old_states = socket_old_states.getDataHandle();
  //activate sub system status to access time, CFL, ...
  Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
  //call function to set time step
  tau=setTimeStep(tau);
  //reset max eigenvalue for computation on next time level
  getMethodData().setMaxEigenval(0.0);

  //creation of copy of actual state for computation of residual after time step
  DataHandle < Framework::State*, Framework::GLOBAL >  states = socket_states.getDataHandle();
  for (CFuint iState = 0; iState < states.size(); ++iState)
  {
    for(CFuint i=0; i<nbEqs; i++)
    {
    (*old_states[iState])[i] = states[iState][0][i];
    }
  }

  //set of time step - Physical time is then updated automatically
  subSysStatus->setDTDim(tau);
  //get matrix of linear systemsolver
  SafePtr<LinearSystemSolver> lss = getMethodData().getLinearSystemSolver()[0];
  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();

  //set trs pointer to inner cells trs
  SafePtr<TopologicalRegionSet> trs = m_cells;
  //get number of inner cells in region inner cells
  const CFuint nbGeos = trs->getLocalNbGeoEnts();

  //get geobuilder to build cells with needed properties (connection, ..)
  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();


  //get structure of data from geobuilder (cells with properties)
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  //set that we use only data of inner cells trs
  geoData.trs = trs;

  //loop over all inner cells
  for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt) {
    CFLogDebugMax("Cell " << iGeoEnt << "\n");

    // build the GeometricEntity (cell)
    //set index of cell
    geoData.idx = iGeoEnt;
    //geo builder make cell
    GeometricEntity& cell = *geoBuilder->buildGE();
    //read vector of states on element (cell)
    std::vector<State*>& cellStates = *(cell.getStates());
    //read number of states on element (cell)
    const CFuint nbStatesInCell = cellStates.size();
    //read vector of nodes on element (cell)

    if (cellStates[0]->isParUpdatable())
    {
      //*****************************************************************
      //*****************************************************************
      //SET Integrator
      //compute shape function in quadrature points
      const std::vector<RealVector>& shapeFunctions =  getMethodData().getVolumeIntegrator()->getSolutionIntegrator(&cell)->computeShapeFunctionsAtQuadraturePoints();

      //numbers of quadrature points
      m_nbKvadrPoint = getMethodData().getVolumeIntegrator()->getSolutionIntegrator(&cell)->getIntegratorPattern()[0];

      //set weights for element quadrature
      const std::valarray<CFreal>& weight = getMethodData().getVolumeIntegrator()->getSolutionIntegrator(&cell)->getCoeff();

      //get coordinates of quadrature points
      const std::vector<RealVector>& coord = getMethodData().getVolumeIntegrator()->getSolutionIntegrator(&cell)->getQuadraturePointsCoordinates();

      //compute gradient of shape functions in quadrature points
      std::vector<RealMatrix> gradient = cell.computeSolutionShapeFunctionGradientsInMappedCoordinates(coord);

      //*****************************************************************
      //*****************************************************************

      /// @todo this must be improved. Finding in a map
      /// for every cell is efficiently speaking not acceptable.
      //make block acumulator with size of sum of states in cell
      DGElemTypeData elemData = m_mapElemData.find(nbStatesInCell);
      BlockAccumulator& acc = *elemData.first;
      RealMatrix& elemMat = *elemData.second;
      RealVector& elemVec = *elemData.third;

      elemVec = 0.0;
      //set matrix in blockaccumulator to 0
      elemMat=0.0;
      acc.setValuesM(elemMat);

      //computation of the Jacobi determinant of mapping from refference element to cell
      /// @todo this must be generalized
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

      //loop over kvadrature point on the cell
      for(CFuint kvadrature_point = 0; kvadrature_point < m_nbKvadrPoint; kvadrature_point++ )
      {
        //set elemMat to 0 if isn't
        if (kvadrature_point!=0) elemMat=0.0;

        //computation of state in point of quadrature
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) //loop over members of state - set to zero
        {
          (*m_state)[iEq] = 0.;
        }
        for (CFuint iState = 0; iState < nbStatesInCell; ++iState) //loop over states in cell
        {
          RealVector &states = *cellStates[iState]->getData();
          for (CFuint iEq = 0; iEq < nbEqs; ++iEq) //loop over members of state
          {
            (*m_state)[iEq] += shapeFunctions[kvadrature_point][iState]*states[iEq];
          }
        }

        //add local rhs to global rhs
        for (CFuint iState = 0; iState < nbStatesInCell; ++iState)
        {
          const CFuint stateID = cellStates[iState]->getLocalID();
          for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
          {
            rhs(stateID, iEq, nbEqs) += (*m_state)[iEq]*shapeFunctions[kvadrature_point][iState]*weight[kvadrature_point]/tau*detJacobi;
          }
        }

        //computation of A_matrix of the cell in point of kvadrature
        if (nbDim == 2)
        {
          compute_Amatrix2D((*m_state),&m_aMatrix);
        }
        else
        {
          compute_Amatrix3D((*m_state),&m_aMatrix);
        }



        for(CFuint s = 0; s < nbDim; s++ ) //loop over index 's' of a matrices
        {
          //compute inner face term
          //loop over test function
          for(CFuint row = 0; row < nbStatesInCell; row++ )
          {
            //loop over base function of solution
            for(CFuint col = 0; col < nbStatesInCell; col++ )
            {
              for(CFuint i = 0; i < nbEqs; i++ )
                for(CFuint j = 0; j < nbEqs; j++ )
                {
                  //inviscid part
                  elemMat(row*nbEqs + i, col*nbEqs + j) -= (m_aMatrix[s](i,j))*gradient[kvadrature_point](row,s)*shapeFunctions[kvadrature_point][col];
                }
            }
          }
        }
        //loop over test function
        for(CFuint row = 0; row < nbStatesInCell; row++ )
        {
          //loop over base function of solution
          for(CFuint col = 0; col < nbStatesInCell; col++ )
          {
            for(CFuint i = 0; i < nbEqs; i++ )
            {
              elemMat(row*nbEqs + i, col*nbEqs + i)+=shapeFunctions[kvadrature_point][row]*shapeFunctions[kvadrature_point][col]/tau;
            }
          }
        }

        //finaly multiply by quadrature weight
        elemMat*=weight[kvadrature_point]*detJacobi;
        // add local matrix to matrix of linear solver using block accumulator
        acc.addValuesM(elemMat);
      }
  //       acc.printToScreen();
      // add the values in the jacobian matrix
      jacobMatrix->addValues(acc);
    }
    //release the GeometricEntity
    geoBuilder->releaseGE();
  }
  CFout << " ... OK\n" << CFendl;

//  jacobMatrix->finalAssembly();
//  jacobMatrix->printToFile("inside");
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdSolveCells::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_states);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > >
StdSolveCells::providesSockets()
{
  std::vector< Common::SafePtr< BaseDataSocketSource > > result;
  result.push_back(&socket_old_states);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

CFreal StdSolveCells::setTimeStep(CFreal tau)
{
  CFAUTOTRACE;
  CFreal norm_square=0.0;
  CFuint m_nbKvadrPoint;

  // get rhs
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  //set trs pointer to inner cells trs
  SafePtr<TopologicalRegionSet> trs = m_cells;
  DataHandle<State*> old_states = socket_old_states.getDataHandle();
  //get number of equations
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  //get dimension
  const CFuint nbDim = PhysicalModelStack::getActive()->getNbEq();
  //get number of inner cells in region inner cells
  const CFuint nbGeos = trs->getLocalNbGeoEnts();

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
    CFLogDebugMax("Cell " << iGeoEnt << "\n");

    // build the GeometricEntity (cell)
    //set index of cell
    geoData.idx = iGeoEnt;
    //geo builder make cell
    GeometricEntity& cell = *geoBuilder->buildGE();

//*****************************************************************
//*****************************************************************
    //SET Integrator
    //compute shape function in quadrature points
    const std::vector<RealVector>& shapeFunctions =  getMethodData().getVolumeIntegrator()->getSolutionIntegrator(&cell)->computeShapeFunctionsAtQuadraturePoints();

    //numbers of quadrature points
    m_nbKvadrPoint = getMethodData().getVolumeIntegrator()->getSolutionIntegrator(&cell)->getIntegratorPattern()[0];

    //get weights for element quadrature
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
//     std::vector<Node*>&  nodes  = *cell.getNodes();


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

    //loop over kvadrature point on the cell
    for(CFuint kvadrature_point = 0; kvadrature_point < m_nbKvadrPoint; kvadrature_point++ )
    {

      for (CFuint iEq = 0; iEq < nbEqs; ++iEq) //loop over members of state - set to zero
      {

        (*m_state)[iEq] = 0;
        (*m_oldState)[iEq] = 0;
      }
      for (CFuint iState = 0; iState < nbStatesInCell; ++iState) //loop over states in cell
      {

        RealVector &state = *cellStates[iState]->getData();
        const CFuint stateID = cellStates[iState]->getLocalID();
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) //loop over members of state
        {

          (*m_state)[iEq] += shapeFunctions[kvadrature_point][iState]*state[iEq];
          (*m_oldState)[iEq] += shapeFunctions[kvadrature_point][iState]*(*old_states[stateID])[iEq];
        }
      }

      //add square of local L2 norm to square of global L2 norm
      for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
      {
        norm_square += ((*m_state)[iEq] - (*m_oldState)[iEq])*((*m_state)[iEq] - (*m_oldState)[iEq])*weight[kvadrature_point]*detJacobi;
      }
    }
    //release the GeometricEntity
    geoBuilder->releaseGE();
  }
  //store norm_1/tau_1 for further use
  CFreal new_Tau;
  static CFreal firstCFL= getMethodData().getCFL()->getCFLValue();
  if (tau == 0)
  {
    new_Tau = firstCFL/(getMethodData().getMaxEigenval());
  }
  else
  {
    static CFreal first_norm = sqrt(norm_square)/tau;
    CFreal newCFL = firstCFL*pow(first_norm/(sqrt(norm_square)/tau),m_Alpha);
    if (m_MaxCFL < newCFL)
    {
      newCFL = m_MaxCFL;
    }
    new_Tau = newCFL/(getMethodData().getMaxEigenval());
    getMethodData().getCFL()->setCFLValue(newCFL);
//     CFout << "     Residual: " << sqrt(norm_square)/tau << "    " << CFendl;
    getMethodData().setResidual(sqrt(norm_square)/tau);
  }
//   Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
// CFout << "\n***** " << subSysStatus->getAllResiduals() << "\n" << CFendl;
// subSysStatus->setResidual(sqrt(norm_square));
// CFout << "\n***** " << subSysStatus->getAllResiduals() << "\n" << CFendl;
//  CFout << sqrt(norm_square) <<  " " << getMethodData().getMaxEigenval() <<  " " << "\n" << CFendl;
//   Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
//   subSysStatus->setAllResiduals(sqrt(norm_square)/tau);
  return(new_Tau);
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

