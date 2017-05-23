#include "Framework/BadFormatException.hh"
#include "Framework/DomainModel.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "FluxReconstructionMethod/FluxReconstruction.hh"
#include "FluxReconstructionMethod/BCStateComputer.hh"

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
  m_faceID(),
  m_trsNames(),
  m_extraVars(CFNULL),
  m_useDomainModel()
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

void BCStateComputer::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  // configure this object by calling the parent class configure()
  FluxReconstructionSolverStrategy::configure(args);
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
//   mapCoordFaceVertexNode[0][KSI] = 0.5;
//   mapCoordFaceVertexNode[1][KSI] = 0.5;

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
          // get boundary face ID
//           const CFuint faceID = face->getID();

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

