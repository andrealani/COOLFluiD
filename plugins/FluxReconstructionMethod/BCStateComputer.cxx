#include "Common/NotImplementedException.hh"
#include "Framework/DiffusiveVarSet.hh"
#include "Framework/BadFormatException.hh"
#include "Framework/DomainModel.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "MathTools/MathFunctions.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "Framework/GeometricEntity.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"
#include "FluxReconstructionMethod/FluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

void BCStateComputer::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("UseDomainModel","Boolean telling whether to use the domain model for the computation of the normals.");
}

//////////////////////////////////////////////////////////////////////////////

BCStateComputer::BCStateComputer(const std::string& name) :
  FluxReconstructionSolverStrategy(name),
  m_needsSpatCoord(),
  m_needsExtraVars(false),
  m_trsNames(),
  m_extraVars(CFNULL),
  m_useDomainModel(),
  m_face(CFNULL),
  m_transitionCriterion(),
  m_nbrTransitionFlags()
{
  CFAUTOTRACE;

  addConfigOptionsTo(this);

  m_useDomainModel = false;
  setParameter("UseDomainModel",&m_useDomainModel);
}

//////////////////////////////////////////////////////////////////////////////

BCStateComputer::~BCStateComputer()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void BCStateComputer::computeBndGradVars(const std::vector< RealVector* >& gradVarsFlxPnt,
                                         const std::vector< Framework::State* >& intStates,
                                         const std::vector< Framework::State* >& ghostStates,
                                         const std::vector< RealVector >& unitNormals,
                                         const std::vector< RealVector >& flxPntCoords,
                                         std::vector< RealVector* >& bndGradVars)
{
  const CFuint nbrStates = intStates.size();
  cf_assert(nbrStates <= gradVarsFlxPnt.size());
  cf_assert(nbrStates <= bndGradVars.size());

  if (nbrStates == 0)
  {
    return;
  }

  const CFuint nbrGradVars = gradVarsFlxPnt[0]->size();

  if (m_gradVarsFace.nbRows() != nbrGradVars || m_gradVarsFace.nbCols() < nbrStates)
  {
    m_gradVarsFace.resize(nbrGradVars,nbrStates);
    m_gradVarsGhost.resize(nbrGradVars,nbrStates);
  }
  m_gradVarStatePtrs.resize(nbrStates);

  SafePtr< DiffusiveVarSet > diffVarSet = getMethodData().getDiffusiveVar();

  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    m_gradVarStatePtrs[iState] = intStates[iState]->getData();
  }
  diffVarSet->setGradientVars(m_gradVarStatePtrs,m_gradVarsFace,nbrStates);

  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    m_gradVarStatePtrs[iState] = ghostStates[iState]->getData();
  }
  diffVarSet->setGradientVars(m_gradVarStatePtrs,m_gradVarsGhost,nbrStates);

  // g_b = a + 0.5*(g(U_ghost) - g(U))
  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    const RealVector& gradVars = *gradVarsFlxPnt[iState];
    RealVector& bndGradVarsFlxPnt = *bndGradVars[iState];

    for (CFuint iVar = 0; iVar < nbrGradVars; ++iVar)
    {
      bndGradVarsFlxPnt[iVar] = gradVars[iVar] + 0.5*(m_gradVarsGhost(iVar,iState) - m_gradVarsFace(iVar,iState));
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCStateComputer::computeBndStates(const std::vector< Framework::State* >& intStates,
                                       const std::vector< Framework::State* >& ghostStates,
                                       const std::vector< RealVector >& unitNormals,
                                       const std::vector< RealVector >& flxPntCoords,
                                       std::vector< RealVector* >& bndStates)
{
  const CFuint nbrStates = intStates.size();

  for (CFuint iState = 0; iState < nbrStates; ++iState)
  {
    const CFuint sizeState = intStates[iState]->size();

    for (CFuint iEq = 0; iEq < sizeState; ++iEq)
    {
      (*bndStates[iState])[iEq] = 0.5*((*(intStates[iState]))[iEq] + (*(ghostStates[iState]))[iEq]);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCStateComputer::computeBndGrads(const std::vector< std::vector< RealVector* > >& intGrads,
                                      std::vector< std::vector< RealVector* > >& bndGrads,
                                      const std::vector< RealVector* >& bndStates,
                                      const std::vector< RealVector >& unitNormals,
                                      const std::vector< RealVector >& flxPntCoords)
{
  // the ghost gradients of this boundary condition, written into bndGrads
  computeGhostGradients(intGrads,bndGrads,unitNormals,flxPntCoords);

  const CFuint nbrFlxPnts = intGrads.size();
  cf_assert(nbrFlxPnts == bndGrads.size());

  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    const CFuint nbrGradVars = intGrads[iFlx].size();

    for (CFuint iVar = 0; iVar < nbrGradVars; ++iVar)
    {
      RealVector& grad = *bndGrads[iFlx][iVar];
      grad = 0.5*(*intGrads[iFlx][iVar] + grad);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCStateComputer::copyGradients(const std::vector< std::vector< RealVector* > >& intGrads,
                                    std::vector< std::vector< RealVector* > >& bndGrads)
{
  const CFuint nbrFlxPnts = intGrads.size();

  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    const CFuint nbrGradVars = intGrads[iFlx].size();

    for (CFuint iVar = 0; iVar < nbrGradVars; ++iVar)
    {
      *bndGrads[iFlx][iVar] = *intGrads[iFlx][iVar];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCStateComputer::removeNormalComponent(RealVector& grad, const RealVector& normal)
{
  const CFreal normalGrad = MathTools::MathFunctions::innerProd(grad,normal);
  grad -= normalGrad*normal;
}

//////////////////////////////////////////////////////////////////////////////

void BCStateComputer::setSlipWallBndGradVars(const std::vector< RealVector* >& gradVarsFlxPnt,
                                             const std::vector< RealVector >& unitNormals,
                                             const std::vector< CFuint >& velocityIDs,
                                             std::vector< RealVector* >& bndGradVars)
{
  const CFuint nbrFlxPnts = gradVarsFlxPnt.size();
  const CFuint nbrVelocities = velocityIDs.size();

  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    *bndGradVars[iFlx] = *gradVarsFlxPnt[iFlx];

    // normal velocity u.n
    CFreal normalVel = 0.;
    for (CFuint iDim = 0; iDim < nbrVelocities; ++iDim)
    {
      normalVel += (*gradVarsFlxPnt[iFlx])[velocityIDs[iDim]]*unitNormals[iFlx][iDim];
    }

    for (CFuint iDim = 0; iDim < nbrVelocities; ++iDim)
    {
      (*bndGradVars[iFlx])[velocityIDs[iDim]] -= normalVel*unitNormals[iFlx][iDim];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCStateComputer::setSlipWallBndGrads(const std::vector< std::vector< RealVector* > >& intGrads,
                                          std::vector< std::vector< RealVector* > >& bndGrads,
                                          const std::vector< RealVector >& unitNormals,
                                          const std::vector< CFuint >& velocityIDs)
{
  copyGradients(intGrads,bndGrads);

  const CFuint nbrFlxPnts = intGrads.size();
  const CFuint nbrVelocities = velocityIDs.size();

  for (CFuint iFlx = 0; iFlx < nbrFlxPnts; ++iFlx)
  {
    const RealVector& normal = unitNormals[iFlx];

    // normal derivative of the normal velocity, n.J.n
    CFreal normalGradUn = 0.;
    for (CFuint iDim = 0; iDim < nbrVelocities; ++iDim)
    {
      normalGradUn += normal[iDim]*MathTools::MathFunctions::innerProd(*intGrads[iFlx][velocityIDs[iDim]],normal);
    }

    // every gradient loses its normal component
    const CFuint nbrGradVars = intGrads[iFlx].size();
    for (CFuint iVar = 0; iVar < nbrGradVars; ++iVar)
    {
      removeNormalComponent(*bndGrads[iFlx][iVar],normal);
    }

    // the velocity block keeps (n.J.n) n n^T
    for (CFuint iDim = 0; iDim < nbrVelocities; ++iDim)
    {
      *bndGrads[iFlx][velocityIDs[iDim]] += (normal[iDim]*normalGradUn)*normal;
    }

    // and loses n (t.grad(u.n)), the tangential gradient of the normal velocity
    m_tangentialGradUn = 0.;
    for (CFuint iDim = 0; iDim < nbrVelocities; ++iDim)
    {
      m_tangentialGradUn += normal[iDim]*(*intGrads[iFlx][velocityIDs[iDim]]);
    }
    removeNormalComponent(m_tangentialGradUn,normal);

    for (CFuint iDim = 0; iDim < nbrVelocities; ++iDim)
    {
      *bndGrads[iFlx][velocityIDs[iDim]] -= normal[iDim]*m_tangentialGradUn;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void BCStateComputer::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  FluxReconstructionSolverStrategy::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void BCStateComputer::setTransitionCriterion(const CFuint iFlux, const bool transition)
{
  cf_assert(m_face != CFNULL);
  std::vector< bool >& flags = m_transitionCriterion[m_face->getID()];
  if (flags.size() != m_nbrTransitionFlags)
  {
    flags.assign(m_nbrTransitionFlags,false);
  }
  flags[iFlux] = transition;
}

//////////////////////////////////////////////////////////////////////////////

bool BCStateComputer::transitionCriterion(const CFuint iFlux) const
{
  cf_assert(m_face != CFNULL);
  std::map< CFuint, std::vector< bool > >::const_iterator flags = m_transitionCriterion.find(m_face->getID());
  return (flags == m_transitionCriterion.end()) ? false : flags->second[iFlux];
}

//////////////////////////////////////////////////////////////////////////////

void BCStateComputer::setup()
{
  CFAUTOTRACE;

  // call setup of parent class
  FluxReconstructionSolverStrategy::setup();

  // add curvature to boundary faces if requested
  /// @note KVDA: this may not be the best place to add curvature to the faces. Consider moving it to StdSetup.
  if (m_useDomainModel)
  {
    addCurvatureToBndFaces();
  }
  
  // get the local FR data
  vector< FluxReconstructionElementData* >& frLocalData = getMethodData().getFRLocalData();
  cf_assert(frLocalData.size() > 0);
  // for now, there should be only one type of element
  cf_assert(frLocalData.size() == 1);
  
  // compute flux point coordinates
  SafePtr< vector<RealVector> > flxLocalCoords = frLocalData[0]->getFaceFlxPntsFaceLocalCoords();
  const CFuint nbrFaceFlxPnts = flxLocalCoords->size();
  
  m_nbrTransitionFlags = nbrFaceFlxPnts;

  // scratch of setSlipWallBndGrads
  m_tangentialGradUn.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

void BCStateComputer::preProcess()
{
}

//////////////////////////////////////////////////////////////////////////////

void BCStateComputer::unsetup()
{
  CFAUTOTRACE;

  // call setup of parent class
  FluxReconstructionSolverStrategy::unsetup();
  
}

//////////////////////////////////////////////////////////////////////////////

void BCStateComputer::addCurvatureToBndFaces()
{
  CFLog(NOTICE, "Adding curvature to the BCs.\n"); 
  if (PhysicalModelStack::getActive()->getDim() == DIM_3D)
  {
    throw Common::NotImplementedException(FromHere(),"Adding of curvature for boundary faces is not yet implemented for 3D...");
  }

  // get the domain model
  SafePtr< DomainModel > domModel = MeshDataStack::getActive()->getDomainModel();

  // get TRS list
  vector< SafePtr< TopologicalRegionSet > > trsList = MeshDataStack::getActive()->getTrsList();

  // number of boundary TRSs with this BC
  const CFuint nbrBCTRSs = m_trsNames.size();

  // get boundary TRSs
  vector< SafePtr< TopologicalRegionSet > > bcTRSs(nbrBCTRSs);
  const CFuint nbTRSs = trsList.size();
  for (CFuint iTRS = 0; iTRS < nbTRSs; ++iTRS)
  {
    for (CFuint iBCTRS = 0; iBCTRS < nbrBCTRSs; ++iBCTRS)
    {
      if (m_trsNames[iBCTRS] == trsList[iTRS]->getName())
      {
        if (bcTRSs[iBCTRS].isNull())
        {
          bcTRSs[iBCTRS] = trsList[iTRS];
        }
        else
        {
          throw BadFormatException (FromHere(),"Two TRSs with the same name found!");
        }
      }
    }
  }

  // some auxiliary variables
  vector< RealVector > mapCoordFaceVertexNode(2,RealVector(1));

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get face builder
  SafePtr< GeometricEntityPool< FaceToCellGEBuilder > > faceBuilder = getMethodData().getFaceBuilder();

  // get the geodata of the face builder and set the cells TRS
  FaceToCellGEBuilder::GeoData& geoData = faceBuilder->getDataGE();
  geoData.cellsTRS = cellTrs;
  geoData.isBoundary = true;

  // loop over boundary condition TRSs
  for (CFuint iBCTRS = 0; iBCTRS < nbrBCTRSs; ++iBCTRS)
  {
    if (bcTRSs[iBCTRS].isNull())
    {
      throw BadFormatException (FromHere(),"Not all boundary TRSs found!");
    }

    // set face TRS in the geodata
    geoData.facesTRS = bcTRSs[iBCTRS];

    // get number of TRs in this TRS
    const CFuint nbTRs = bcTRSs[iBCTRS]->getNbTRs();

    // loop over the TRs
    for (DomainModel::TRidx iTR = 0; iTR < nbTRs; ++iTR)
    {
      // get the TR
      SafePtr< TopologicalRegion > bcTR = bcTRSs[iBCTRS]->getTopologicalRegion(iTR);

      // get TR global index
      const std::string trKey = m_trsNames[iBCTRS] + StringOps::to_str(iTR);
      DomainModel::TRidx trGlobalIdx = domModel->getTRGlobalIdx(trKey);

      // loop over faces in this TR
      const CFuint nbrFaces = bcTR->getLocalNbGeoEnts();
      for (CFuint iFace = 0; iFace < nbrFaces; ++iFace)
      {
        // get face index in the TRS
        const CFuint faceIdx = bcTR->getGeoIDInTrs(iFace);

        // build the face GeometricEntity
        geoData.idx = faceIdx;
        GeometricEntity* face = faceBuilder->buildGE();

        // check if face is higher order than P1
        if (face->getGeometryShapeFunctionOrder() > CFPolyOrder::ORDER1)
        {
          // get face nodes
          vector< Node* >* nodes = face->getNodes();
          cf_assert(nodes->size() > 2);

          // get the mapped coordinates of the face vertex nodes
          mapCoordFaceVertexNode[0] = 0.5; // copy for a good guess
          domModel->computeParamCoord(trGlobalIdx,*(*nodes)[0],mapCoordFaceVertexNode[0]);
          mapCoordFaceVertexNode[1] = 0.5; // copy for a good guess
          domModel->computeParamCoord(trGlobalIdx,*(*nodes)[1],mapCoordFaceVertexNode[1]);

          /// @warning KVDA: the following is hard coded for P2 geometrical order
          // compute mapped coordinate of the `inner' face node
          const RealVector mapCoordFaceInnderNode = 0.5*(mapCoordFaceVertexNode[0]+
              mapCoordFaceVertexNode[1]);

          // compute the new coordinates from the mapped coordinates
          domModel->computeCoord(trGlobalIdx,mapCoordFaceInnderNode,*(*nodes)[2]);

          // if the neighbouring cell is CFGeoShape::QUAD, also move the internal cell node
          GeometricEntity* cell = face->getNeighborGeo(0);
          if (cell->getShape() == CFGeoShape::QUAD)
          {
            // get the cell nodes
            vector< Node* >* nodes = cell->getNodes();
            cf_assert(nodes->size() == 9);

            // recompute the inner cell node from the other nodes
            *(*nodes)[8] = - 0.25*(*(*nodes)[0] + *(*nodes)[1] + *(*nodes)[2] + *(*nodes)[3])
                + 0.50*(*(*nodes)[4] + *(*nodes)[5] + *(*nodes)[6] + *(*nodes)[7]);
          }

          /// @todo add check for Jacobian determinant positivity
//           cell->computeGeometricShapeFunctionJacobianDeterminant(const std::vector<RealVector>& mappedCoord);
        }

        // release the face
        faceBuilder->releaseGE();
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace FluxReconstructionMethod

}  // namespace COOLFluiD

