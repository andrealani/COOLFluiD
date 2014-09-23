#include "MeshLaplacianSmoothing/MeshLaplacianSmoothing.hh"


#include "UpdateMesh.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "MeshTools/QualityCalculator.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshLaplacianSmoothing {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<UpdateMesh, LaplacianSmoothingData, MeshLaplacianSmoothingModule> UpdateMeshProvider("UpdateMesh");

//////////////////////////////////////////////////////////////////////////////

void UpdateMesh::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("NbSmoothingIter","Number of Smoothing Iterations");
   options.addConfigOption< std::string >("Weight","Type of weightfunction to use Uniform/overLength/overDistance");
   options.addConfigOption< CFreal >("Relaxation","Relaxation Factor");
   options.addConfigOption< CFreal >("QualityThreshold","Threshold for selecting which nodes will be moved.");
}

//////////////////////////////////////////////////////////////////////////////

UpdateMesh::UpdateMesh(const std::string& name) :
LaplacianSmoothingCom(name),
  socket_nodes("nodes"),
  socket_averageVector("averageVector"),
  socket_sumWeight("sumWeight"),
  socket_wallDistance("wallDistance"),
  socket_qualityNode("qualityNode"),
  socket_nodeToCellConnectivity("nodeToCellConnectivity"),
  _coordI()
{
   addConfigOptionsTo(this);
  _relaxation = 0.5;
   setParameter("Relaxation",&_relaxation);

  _nbSmoothingIter = 10;
   setParameter("NbSmoothingIter",&_nbSmoothingIter);

  _weightType = "Uniform";
   setParameter("Weight",&_weightType);

  _qualityThreshold = 1.;
   setParameter("QualityThreshold",&_qualityThreshold);

}

//////////////////////////////////////////////////////////////////////////////

void UpdateMesh::setup()
{
  LaplacianSmoothingCom::setup();

  _coordI.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
UpdateMesh::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_averageVector);
  result.push_back(&socket_sumWeight);
  result.push_back(&socket_wallDistance);
  result.push_back(&socket_qualityNode);
  result.push_back(&socket_nodeToCellConnectivity);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void UpdateMesh::execute()
{
  CFAUTOTRACE;
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle< RealVector> averageVector = socket_averageVector.getDataHandle();
  DataHandle< CFreal> sumWeight = socket_sumWeight.getDataHandle();
  DataHandle< CFreal> qualityNode = socket_qualityNode.getDataHandle();
  DataHandle< vector<CFuint> > nodeToCellConn = socket_nodeToCellConnectivity.getDataHandle();

  Common::SelfRegistPtr<MeshTools::QualityCalculator> qualityComputer = Environment::Factory<MeshTools::QualityCalculator>::getInstance().
      getProvider("Concrete")->create("Concrete");
  cf_assert(qualityComputer.isNotNull());

  Common::SafePtr<TopologicalRegionSet> trs =
    MeshDataStack::getActive()->getTrs("InnerCells");

  const CFreal oldQualityThreshold = _qualityThreshold;
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbNodes = nodes.size();
  const CFuint maxNbIter = 1;
  const CFreal oldRelaxation = _relaxation;
  RealVector oldCoord(nbDim);
  bool isAchieved = false;
  bool onlyAllowImprovements = false;

//  std::cout << "Moving the nodes using the LaplacianSmoothing" << std::endl;
  for (_iIter = 0; (!isAchieved); ++_iIter)
  {
//    std::cout << "Smoothing Step: " << iIter << std::endl;
    // Compute the average 'displacement Vector'
    computeDisplacementVectors();

    // Move the nodes
    CFuint nbModifiedNodes = 0;
    CFuint nbNegativeVolume = 0;
    for (CFuint iNode = 0; iNode < nbNodes; ++iNode)
    {
      bool isLocallyAchieved = false;
      CFuint localID = nodes[iNode]->getLocalID();
      for (CFuint iLocalIter = 0; (!isLocallyAchieved); ++iLocalIter)
      {
        if(((qualityNode[localID] > _qualityThreshold)) && (qualityNode[localID] != 0.))
        {
          //Save oldCoord
          oldCoord = (*nodes[iNode]);

          //Save oldQuality
          CFreal oldQuality = qualityNode[localID];

          // Compute the new coordinates
          computeNewCoordinates(iNode);

          //Move Node
          for (CFuint iDim = 0; iDim < nbDim; ++iDim)
          {
            (*nodes[iNode])[iDim] = _coordI[iDim];
          }

          //Check that the worse quality has been improved
          vector<CFuint>& cellsID = nodeToCellConn[localID];

          Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
          geoBuilder = getMethodData().getStdTrsGeoBuilder();

          StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
          geoData.trs = trs;

          CFreal worseQuality = 0.;
          for (CFuint iCell=0; iCell < cellsID.size(); ++iCell) {
            // build the GeometricEntity
            geoData.idx = cellsID[iCell];
            GeometricEntity* currCell = geoBuilder->buildGE();

            CFreal quality = qualityComputer->computeQuality(currCell);
            if (quality > worseQuality) worseQuality = quality;
            if (quality < 0.) worseQuality = MathTools::MathConsts::CFrealMax();

            geoBuilder->releaseGE();
          }

          if((worseQuality > 10.E98)&&(oldQuality > 10.E98))
          {
            _hasNegativeVolumeCells = true;
            nbNegativeVolume++;
          }

          if(((worseQuality > 10.E98)&&(oldQuality < 10.E98))||(worseQuality > oldQuality) && (onlyAllowImprovements))
          {
            //std::cout<< "Trying to create a negative volume cell...Reverting" << std::endl;
            _relaxation *= 0.5;
            (*nodes[iNode]) = oldCoord;
          }
          else
          {
          //std::cout<< "Quality improved from: " << oldQuality << " to " << worseQuality << std::endl;
            qualityNode[localID] = worseQuality;
            _relaxation = oldRelaxation;
            nbModifiedNodes++;
          }
        }
        if(iLocalIter >= maxNbIter) isLocallyAchieved = true;
      }
      averageVector[localID] = 0.;
      sumWeight[localID] = 0.;
    }
    std::cout << "Number of modifications: " << nbModifiedNodes << std::endl;
    if(_iIter >= _nbSmoothingIter) isAchieved = true;

    if(nbNegativeVolume > 0){
      _hasNegativeVolumeCells = true;
      std::cout<< "We still have " << nbNegativeVolume << " negative volumes" << std::endl;
      if ( !(_iIter % 3) ) {
        std::cout<< "Threshold modified from : " <<  _qualityThreshold;
        _qualityThreshold = 1. + ((_qualityThreshold-1.)*0.5);
        std::cout<< " to: " <<  _qualityThreshold << std::endl;
      }

      if(_iIter > _nbSmoothingIter){
        isAchieved = false;
        _relaxation *= 0.5;
        _qualityThreshold = 1.;
      }
    }

    if((_iIter >= 2*_nbSmoothingIter) || ((_hasNegativeVolumeCells == false)&&(_iIter > _nbSmoothingIter))) isAchieved = true;
  }

  checkNewCells();
  if(_hasNegativeVolumeCells == true)
      cout << "WARNING: The new mesh contains " << _nbNegativeVolumeCells << " cells with a negative volume" << endl;
//  cf_assert(_hasNegativeVolumeCells == false);

  //Set back the quality threshold
  _qualityThreshold = oldQualityThreshold;

}

//////////////////////////////////////////////////////////////////////////////

void UpdateMesh::computeNewCoordinates(const CFuint iNode)
{

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle< RealVector> averageVector = socket_averageVector.getDataHandle();
  DataHandle< CFreal> sumWeight = socket_sumWeight.getDataHandle();

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint localID = (*nodes[iNode]).getLocalID();
  for (CFuint iDim = 0; iDim < nbDim; ++iDim)
  {
    _coordI[iDim] = (1.-_relaxation) * (*nodes[iNode])[iDim];
    _coordI[iDim] += _relaxation * (averageVector[localID][iDim]/sumWeight[localID]);
  }

}

//////////////////////////////////////////////////////////////////////////////

void UpdateMesh::computeDisplacementVectors()
{
  DataHandle< RealVector> averageVector = socket_averageVector.getDataHandle();
  DataHandle< CFreal> sumWeight = socket_sumWeight.getDataHandle();
  DataHandle<CFreal> wallDistance = socket_wallDistance.getDataHandle();

  RealVector edgeVector(PhysicalModelStack::getActive()->getDim());
  RealVector jDisplacement(PhysicalModelStack::getActive()->getDim());
  CFreal weight = 1.0;

  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->
   getTrs("InnerCells");

  const CFuint nbCells = cells->getLocalNbGeoEnts();
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();

  //iNITIALIZE
  for (CFuint iNode = 0; iNode < averageVector.size(); ++iNode)
  {
    averageVector[iNode] = 0;
    sumWeight[iNode] = 0;
  }

  Common::SafePtr<GeometricEntityPool<TrsGeoWithNodesBuilder> >
    geoBuilder = getMethodData().getGeoWithNodesBuilder();

  TrsGeoWithNodesBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = cells;

  for (CFuint iCell = 0; iCell < nbCells; ++iCell)
  {
    // build the GeometricEntity
    geoData.idx = iCell;
    GeometricEntity *const currCell = geoBuilder->buildGE();
    std::vector<Node*>* cellNodes = currCell->getNodes();
    const CFuint nbNodesInCell = cellNodes->size();
    for (CFuint iNode = 0; iNode < nbNodesInCell; ++iNode)
      {
        const CFuint localID = (*cellNodes)[iNode]->getLocalID();
        _coordI = *((*cellNodes)[iNode]);

        for (CFuint jNode = 0; jNode < nbNodesInCell; ++jNode)
        {
          if(jNode != iNode)
          {

            if(_weightType == "overLength")
            {
              edgeVector = *((*cellNodes)[jNode]) - _coordI;
              CFreal length = edgeVector[0]*edgeVector[0];
              for(CFuint iDim=1;iDim <nbDim;iDim++){
                length += edgeVector[iDim]*edgeVector[iDim];
                }
              weight = 1./(sqrt(length) + MathTools::MathConsts::CFrealEps());
            }

            if(_weightType == "Uniform")
            {
              weight = 1.;
            }

            if(_weightType == "overDistance")
            {
              const CFuint jLocalID = (*cellNodes)[jNode]->getLocalID();
              CFreal distance = wallDistance[jLocalID];
              weight = 1./(distance+MathTools::MathConsts::CFrealEps());
            }

            edgeVector = *((*cellNodes)[jNode]);
            edgeVector *= weight;

            averageVector[localID] += edgeVector;
            sumWeight[localID] += weight;
          }
        } // end looping over Nodes of Cell != currentNode
      } // end looping over Nodes of Cell

    //release the GeometricEntity
    geoBuilder->releaseGE();
  } // end looping over Cells
}

//////////////////////////////////////////////////////////////////////////////

void UpdateMesh::checkNewCells()
{

  _hasNegativeVolumeCells = false;

  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  _nbNegativeVolumeCells = 0;

  Common::SafePtr<GeometricEntityPool<TrsGeoWithNodesBuilder> >
    geoBuilder = getMethodData().getGeoWithNodesBuilder();

  TrsGeoWithNodesBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = cells;

  // Compute the cell Volume
  for (CFuint iCell = 0; iCell < nbCells; ++iCell)
  {
    // build the GeometricEntity
    geoData.idx = iCell;
    GeometricEntity *const currCell = geoBuilder->buildGE();

    CFreal cellVolume = currCell->computeVolume();
    if(cellVolume < 0.)
    {
     _nbNegativeVolumeCells++;
    }

    //release the GeometricEntity
    geoBuilder->releaseGE();
  }

  if(_nbNegativeVolumeCells > 0){
    _hasNegativeVolumeCells = true;
  }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshLaplacianSmoothing

  } // namespace Numerics

} // namespace COOLFluiD
