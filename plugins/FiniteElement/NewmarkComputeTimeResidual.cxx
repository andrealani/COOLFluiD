#include "Framework/MeshData.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"

#include "FiniteElement/FiniteElement.hh"
#include "FiniteElement/NewmarkComputeTimeResidual.hh"
#include "FiniteElement/ComputeInertiaTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NewmarkComputeTimeResidual, FiniteElementMethodData, FiniteElementModule> newmarkComputeTimeResidualProvider("NewmarkComputeTimeResCom");

//////////////////////////////////////////////////////////////////////////////

void NewmarkComputeTimeResidual::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Alpha","Alpha coef");
   options.addConfigOption< CFreal >("Gamma","Gamma coef");
   options.addConfigOption< bool >("LumpedMassMatrix","Use or not a lumped Lumped Mass Matrix using sum row technique");
}

//////////////////////////////////////////////////////////////////////////////

NewmarkComputeTimeResidual::NewmarkComputeTimeResidual(const std::string& name) :
FiniteElementMethodCom(name),
  socket_rhs("rhs"),
  socket_pastStates("pastStates"),
  socket_pastStatesD("pastStatesD"),
  socket_pastStatesD2("pastStatesD2"),
  socket_appliedStrongBC("appliedStrongBC")
{
   addConfigOptionsTo(this);
  _alpha = 0.5;
   setParameter("Alpha",&_alpha);

  _gamma = 0.5;
   setParameter("Gamma",&_gamma);

  _lumpMassMatrix = false;
   setParameter("LumpedMassMatrix",&_lumpMassMatrix);

}

//////////////////////////////////////////////////////////////////////////////

NewmarkComputeTimeResidual::~NewmarkComputeTimeResidual()
{
}

//////////////////////////////////////////////////////////////////////////////

void NewmarkComputeTimeResidual::setup()
{
  CFAUTOTRACE;

  // first call parent method
  FiniteElementMethodCom::setup();

  // rhs storage
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  SafePtr<vector<ElementTypeData> > elementType =
    MeshDataStack::getActive()->getElementTypeData();

  const CFuint nbElemTypes = elementType->size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  // settin correct size for integration result
   _integResult.resize(nbEqs,nbEqs);

  // loop over types since it can happen to deal with an hybrid mesh
   for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
     // number of states in a cell of this type
     const CFuint nbStatesInType = (*elementType)[iType].getNbStates();
     // create a block matrix accumulator for this type of cell
     BlockAccumulator* ptr = getMethodData().getLinearSystemSolver()[0]->
      createBlockAccumulator(nbStatesInType,nbStatesInType,nbEqs);
    _mapAcc.insert(nbStatesInType,ptr);
   }

  _mapAcc.sortKeys();
}

//////////////////////////////////////////////////////////////////////////////

void NewmarkComputeTimeResidual::executeOnTrs()
{
  CFAUTOTRACE;

  // rhs storage
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  // storage of the past States
  DataHandle<State*> pastStates  = socket_pastStates.getDataHandle();

  // storage of the past States derivatives
  DataHandle<State*> pastStatesD  = socket_pastStatesD.getDataHandle();

  // storage of the past States second derivatives
  DataHandle<State*> pastStatesD2  = socket_pastStatesD2.getDataHandle();

  //
  DataHandle<std::vector<bool> > appliedStrongBC  = socket_appliedStrongBC.getDataHandle();

  SafePtr<LinearSystemSolver> lss = getMethodData().getLinearSystemSolver()[0];
  SafePtr<ComputeInertiaTerm> inertiaComputer = getMethodData().getInertiaTermComputer();

  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();

  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbGeos = trs->getLocalNbGeoEnts();
  FiniteElementMethodData& femdata  = getMethodData();
  LocalElementData& local_elem_data = femdata.getLocalElementData();
  local_elem_data.nbEqs = nbEqs;
  local_elem_data.trs   = getCurrentTRS();
  PhysicalModelStack::getActive()->setCurrentZone(getCurrentTRS()->getName());

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
    geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = trs;

  const CFreal dt = SubSystemStatusStack::getActive()->getDT();
  const CFreal a3 = 2./(_gamma * dt * dt);
  const CFreal a4 = dt * a3;
  const CFreal a5 = (1./_gamma) - 1.;

  std::vector<RealVector> inertiaRhs(0);
  RealMatrix massMatrix(nbEqs,nbEqs);

  SafePtr<BlockAccumulator> acc;

  for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt) {

    CFLogDebugMax("Cell " << iGeoEnt << "\n");

    // build the GeometricEntity
    geoData.idx = iGeoEnt;
    GeometricEntity *const cell = geoBuilder->buildGE();
    vector<State*>* cellStates = cell->getStates();

    getMethodData().getFEMVolumeIntegrator()->precomputeCellData(cell);
//   BlockAccumulator& acc = *local_elem_data.blockacc;
//  RealMatrix& elemMat   = *local_elem_data.stiff_mat;
//   RealVector& elemVec   = *local_elem_data.load_vec;

    local_elem_data.nbStates = cellStates->size();
    local_elem_data.cell = cell;


    /// @todo this must be improved. Finding in a map for every cell is effieciently speaking not acceptable.
      acc = _mapAcc.find(cellStates->size());
      acc->reset();

      // Compute and  Insert result into the Block Accumulator
      CFuint nbStatesInCell = cell->getNbNodesSolutionShapeFunction();

      inertiaRhs.resize(nbStatesInCell);
      for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
        inertiaRhs[iState].resize(nbEqs);
      }

      /// @todo AL: this is hell inefficient !!!
      RealMatrix elemMat(nbStatesInCell*nbEqs,nbStatesInCell*nbEqs);

      for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
        const CFuint stateID = (*cellStates)[iState]->getLocalID();
        acc->setRowColIndex(iState, stateID);
      }

      for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
        State *const currState = (*cellStates)[iState];
        CFuint stateID = currState->getLocalID();
        local_elem_data.iState = iState;

//if one of the variable is dirichlet, then dont put the state!!!
///@todo change this??!!!
bool dirichlet(false);
if(appliedStrongBC[stateID].size() > 0) dirichlet = true;

        // Computing the part of the Inertia term going to the RHS (in the explicit case)
        // This term is: M *(a3*u_n + a4*up_n + a5*upp_n)
        // so we first compute: (a3*u_n + a4*up_n + a5*upp_n)
        inertiaRhs[iState] =  a3 * (*pastStates[stateID]);
        inertiaRhs[iState] += a4 * (*pastStatesD[stateID]);
        inertiaRhs[iState] += a5 * (*pastStatesD2[stateID]);

        if (currState->isParUpdatable() && !dirichlet) {
          for (CFuint jState = 0; jState < nbStatesInCell; ++jState) {

            // Inertia Term
            local_elem_data.jState = jState;
            inertiaComputer->computeTerm(cell, _integResult);

            // Newmark method, multiply by a3
            massMatrix = a3 * _integResult;

            if(_lumpMassMatrix){
              // Adding part of Mass matrix going to system matrix
              // with Row-Summing Lumping
              CFLogDebugMax("Adding values to " << iState << " " << iState << "\n");
              acc->addValuesM(iState,
			      iState,
			      massMatrix);

              // Compute Element Mass Matrix
              // with Row-Summing Lumping
              for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
                for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
                  elemMat(iState*nbEqs + iEq, iState*nbEqs + jEq) += _integResult(iEq,jEq);
                }
              }
            }
            else{
              // Adding part of Mass matrix going to system matrix
              // without Row-Summing Lumping
              CFLogDebugMax("Adding values to " << iState << " " << jState << "\n");
              acc->addValuesM(iState,
			      jState,
			      massMatrix);

              // Compute Element Mass Matrix
              // with Row-Summing Lumping
              for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
                for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
                  elemMat(iState*nbEqs + iEq, jState*nbEqs + jEq) += _integResult(iEq,jEq);
                }
              }
            }

          }
        }
      }

      // write to the right hand side
      /// @todo try to improve this product
      for (CFuint iState = 0; iState < nbStatesInCell; ++iState) {
        CFuint stateID = (*cellStates)[iState]->getLocalID();
        if(appliedStrongBC[stateID].size() == 0)
        {
          for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
            for (CFuint jState = 0; jState < nbStatesInCell; ++jState) {
              RealVector& jStateRhs = inertiaRhs[jState];
              for (CFuint jEq = 0; jEq < nbEqs; ++jEq) {
                rhs(stateID, iEq, nbEqs) += elemMat(iState*nbEqs + iEq, jState*nbEqs + jEq) * jStateRhs[jEq];
              }
            }
          }
        }
      }

      // add the values in the jacobian matrix
      jacobMatrix->addValues(*acc);

      //release the GeometricEntity
      geoBuilder->releaseGE();
  }
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
NewmarkComputeTimeResidual::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);
  result.push_back(&socket_pastStates);
  result.push_back(&socket_pastStatesD);
  result.push_back(&socket_pastStatesD2);
  result.push_back(&socket_appliedStrongBC);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD
