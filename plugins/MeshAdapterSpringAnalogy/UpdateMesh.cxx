#include "MeshAdapterSpringAnalogy/MeshAdapterSpringAnalogy.hh"


#include "UpdateMesh.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshAdapterSpringAnalogy {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<UpdateMesh, SpringAnalogyData, MeshAdapterSpringAnalogyModule> UpdateMeshProvider("UpdateMesh");

//////////////////////////////////////////////////////////////////////////////

void UpdateMesh::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("NbSmoothingIter","Number of Smoothing Iterations");
   options.addConfigOption< std::string >("Weight","Type of weightfunction to use Uniform/Length/overLength/Pirzadeh");
   options.addConfigOption< CFreal >("Relaxation","Relaxation Factor");
}

//////////////////////////////////////////////////////////////////////////////

UpdateMesh::UpdateMesh(const std::string& name) :
SpringAnalogyCom(name),
  socket_nodes("nodes"),
  socket_averageVector("averageVector"),
  socket_sumWeight("sumWeight"),
  socket_isMovable("isMovable"),
  socket_nodalDisplacements("nodalDisplacements"),
  socket_wallDistance("wallDistance",false),
  _coordI()
{
   addConfigOptionsTo(this);
  _relaxation = 1.;
   setParameter("Relaxation",&_relaxation);

  _nbSmoothingIter = 10;
   setParameter("NbSmoothingIter",&_nbSmoothingIter);

  _weightType = "Uniform";
   setParameter("Weight",&_weightType);

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
UpdateMesh::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_averageVector);
  result.push_back(&socket_sumWeight);
  result.push_back(&socket_isMovable);
  result.push_back(&socket_nodalDisplacements);
  result.push_back(&socket_wallDistance);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void UpdateMesh::setup()
{
  SpringAnalogyCom::setup();

  _coordI.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

void UpdateMesh::execute()
{
  CFAUTOTRACE;
  
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle< RealVector> averageVector = socket_averageVector.getDataHandle();
  DataHandle< CFreal> sumWeight = socket_sumWeight.getDataHandle();
  DataHandle< bool> isMovable = socket_isMovable.getDataHandle();
  DataHandle< RealVector> displacements = socket_nodalDisplacements.getDataHandle();

  _ballVertexComputer.setup();

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbNodes = nodes.size();

  RealVector oldCoord(nbDim);
  bool isAchieved = false;

  for (_iIter = 0; (!isAchieved); ++_iIter)
  {
    // Compute the average 'displacement Vector'
    computeDisplacementVectors();

    // Move the nodes
    for (CFuint iNode = 0; iNode < nbNodes; ++iNode)
    {
      CFuint localID = nodes[iNode]->getLocalID();

      if(isMovable[localID] != false)
      {
        //Save oldCoord
        oldCoord = (*nodes[iNode]);

        // Compute the new coordinates
        computeNewCoordinates(iNode);

        //Move Node
        for (CFuint iDim = 0; iDim < nbDim; ++iDim)
        {
          (*nodes[iNode])[iDim] = _coordI[iDim];
        }
        displacements[localID] = (*nodes[iNode]) - oldCoord;
      }

      averageVector[localID] = 0.;
      sumWeight[localID] = 0.;
    }

    // Set the displacement to zero for the boundary nodes
    for (CFuint iNode = 0; iNode < nbNodes; ++iNode)
    {
      CFuint localID = nodes[iNode]->getLocalID();
      if(isMovable[localID] != true)
      {
        displacements[localID] = 0;
      }
    }

    if(_iIter >= _nbSmoothingIter) isAchieved = true;
  }

  CFout << "Checking for negative cells "<< "\n" << CFendl;

  checkNewCells();

  if(_hasNegativeVolumeCells == true)
      cout << "WARNING: The new mesh contains " << _nbNegativeVolumeCells << " cells with a negative volume" << endl;
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
      _coordI[iDim] = (*nodes[iNode])[iDim];
      _coordI[iDim] += _relaxation * (averageVector[localID][iDim]/sumWeight[localID]);
  }

}

//////////////////////////////////////////////////////////////////////////////

void UpdateMesh::computeDisplacementVectors()
{
  RealVector edgeVector(PhysicalModelStack::getActive()->getDim());
  CFreal weight = 0.0;

  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  DataHandle< RealVector> averageVector = socket_averageVector.getDataHandle();
  DataHandle< CFreal> sumWeight = socket_sumWeight.getDataHandle();
  DataHandle< RealVector> displacements = socket_nodalDisplacements.getDataHandle();
  DataHandle< bool> isMovable = socket_isMovable.getDataHandle();

  const CFuint nbCells = cells->getLocalNbGeoEnts();

  Common::SafePtr<GeometricEntityPool<TrsGeoWithNodesBuilder> >
    geoBuilder = getMethodData().getGeoWithNodesBuilder();

  TrsGeoWithNodesBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = cells;

  for (CFuint iNode = 0; iNode < averageVector.size(); ++iNode) {
    averageVector[iNode] = 0;
    sumWeight[iNode] = 0;
  }

  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    // build the GeometricEntity
    geoData.idx = iCell;
    GeometricEntity *const currCell = geoBuilder->buildGE();

    std::vector<Node*>* cellNodes = currCell->getNodes();
    const CFuint nbNodesInCell = cellNodes->size();
    for (CFuint iNode = 0; iNode < nbNodesInCell; ++iNode) {
      const CFuint localID = (*cellNodes)[iNode]->getLocalID();
      _coordI = *((*cellNodes)[iNode]);

      //          CFreal sinA;
      //          if(_weightType == "Murayama") sinA = computeNodalAngle(cellNodes, iNode);

      bool closeBoundary = false;
      for (CFuint jNode = 0; jNode < nbNodesInCell; ++jNode) {
	if(jNode != iNode) {
	  const CFuint jLocalID = (*cellNodes)[jNode]->getLocalID();
	  if(isMovable[jLocalID] == false) closeBoundary = true;
	}
      }

      for (CFuint jNode = 0; jNode < nbNodesInCell; ++jNode) {
	if(jNode != iNode) {
	  const CFuint jLocalID = (*cellNodes)[jNode]->getLocalID();

	  if(_weightType == "overDistance") {
            DataHandle<CFreal> wallDistance = socket_wallDistance.getDataHandle();
	    const CFreal distance = wallDistance[jLocalID];
	    weight = 1./(distance+MathTools::MathConsts::CFrealEps());
	    edgeVector = displacements[jLocalID];
	    edgeVector *= weight;
	  }

	  if(_weightType == "overDistance2") {
            DataHandle<CFreal> wallDistance = socket_wallDistance.getDataHandle();
	    const CFreal distance = wallDistance[jLocalID];
	    weight = 1./(distance*distance+MathTools::MathConsts::CFrealEps());
	    edgeVector = displacements[jLocalID];
	    edgeVector *= weight;
	  }

	  if(_weightType == "overDistance2Link") {
            DataHandle<CFreal> wallDistance = socket_wallDistance.getDataHandle();
	    CFreal distance = wallDistance[jLocalID];
            distance /= 150.;
	    if ((closeBoundary == true && isMovable[jLocalID] != true)
		||(closeBoundary == false)){
	      weight = 1./(distance*distance+MathTools::MathConsts::CFrealEps());
	      edgeVector = displacements[jLocalID];
	      //if(distance < distanceI) weight *= 2.;
	      edgeVector *= weight;
	    }
	  }

	  if(_weightType == "overDistance3") {
            DataHandle<CFreal> wallDistance = socket_wallDistance.getDataHandle();
	    const CFreal distance = wallDistance[jLocalID];
            cf_assert (distance > 0.);
	    weight = 1./(distance*distance*distance+MathTools::MathConsts::CFrealEps());
	    edgeVector = displacements[jLocalID];
	    edgeVector *= weight;
	  }

	  if(_weightType == "overDistance4") {
            DataHandle<CFreal> wallDistance = socket_wallDistance.getDataHandle();
	    const CFreal distance = wallDistance[jLocalID];
            cf_assert (distance > 0.);
	    weight = 1./(distance*distance*distance*distance+MathTools::MathConsts::CFrealEps());
	    edgeVector = displacements[jLocalID];
	    edgeVector *= weight;
	  }

	  if(_weightType == "overLength") {
	    edgeVector = *((*cellNodes)[jNode]) - _coordI;
	    CFreal length = edgeVector.norm2();
	    weight = 1./(sqrt(length) + MathTools::MathConsts::CFrealEps());

	    edgeVector = displacements[jLocalID];
	    edgeVector *= weight;
	  }

	  if(_weightType == "Uniform") {
	    edgeVector = *((*cellNodes)[jNode]) - _coordI;
	    weight = 1.;

	    edgeVector = displacements[jLocalID];
	    edgeVector *= weight;
	  }

	  //             if(_weightType == "Murayama")
	  //             {
	  //               edgeVector = *((*cellNodes)[jNode]) - _coordI;
	  //               CFreal length = edgeVector[0]*edgeVector[0];
	  //               for(CFuint iDim=1;iDim <nbDim;iDim++){
	  //                 length += edgeVector[iDim]*edgeVector[iDim];
	  //                 }
	  //               weight = 1./(length+MathTools::MathConsts::CFrealEps());
	  //               weight += 1./(sinA*sinA);
	  //               weight = min(weight/wallDistance[localID], 500.);
	  //               CFreal distance = wallDistance[localID];
	  //               if(distance > 0.1) weight = 1.;
	  //             }

	  if(_weightType == "overDistance_BallVertex") {
            DataHandle<CFreal> wallDistance = socket_wallDistance.getDataHandle();
	    edgeVector = *((*cellNodes)[jNode]) - _coordI;
	    //CFreal length = edgeVector.norm2();

	    CFreal distance = wallDistance[jLocalID];
	    weight = 1./(distance+MathTools::MathConsts::CFrealEps());

	    //Linear Spring part
	    edgeVector = displacements[jLocalID];
	    edgeVector *= weight;

	    //Ball Vertex part
	    RealVector ballVertex =
	      _ballVertexComputer.computeBallVertex(currCell, iNode, jNode);

	    CFreal weightBall = 1.;
	    ballVertex *= weightBall;
	    weight += weightBall;
	    edgeVector += ballVertex;
	    //std::cout << edgeVector << std::endl;
	  }

          if(_weightType == "Uniform_BallVertex") {
              //Uniform part
              weight = 1.;

              //Linear Spring part
              edgeVector = displacements[jLocalID];

              //Ball Vertex part
              edgeVector += _ballVertexComputer.computeBallVertex(currCell, iNode, jNode);

              //std::cout << edgeVector << std::endl;
              edgeVector *= weight;
              weight += 1.;
           }

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

  Common::SafePtr<TopologicalRegionSet> cells =
    MeshDataStack::getActive()->getTrs("InnerCells");
  const CFuint nbCells = cells->getLocalNbGeoEnts();
  _nbNegativeVolumeCells = 0;

  Common::SafePtr<GeometricEntityPool<TrsGeoWithNodesBuilder> >
    geoBuilder = getMethodData().getGeoWithNodesBuilder();

  TrsGeoWithNodesBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.trs = cells;

  // Compute the cell Volume
  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    // build the GeometricEntity
    geoData.idx = iCell;
    GeometricEntity *const currCell = geoBuilder->buildGE();
    const CFreal cellVolume = currCell->computeVolume();
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

// CFreal UpdateMesh::computeNodalAngle(std::vector<Node*>* cellNodes, CFuint iNode)
// {
//   CFuint node1;
//   CFuint node2;
//   if(iNode == 0){
//     node1 = 1;
//     node2 = 2;
//     }
//   if(iNode == 1){
//     node1 = 2;
//     node2 = 0;
//     }
//   if(iNode == 2){
//     node1 = 0;
//     node2 = 1;
//     }
//
//   const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
//   RealVector edgeVector(nbDim);
//   edgeVector = *((*cellNodes)[node1]) - *((*cellNodes)[node2]);
//   CFreal length = edgeVector[0]*edgeVector[0];
//   for(CFuint iDim=1;iDim <nbDim;iDim++){
//     length += edgeVector[iDim]*edgeVector[iDim];
//     }
//
//   CFreal a = sqrt(length);
//
//   edgeVector = *((*cellNodes)[iNode]) - *((*cellNodes)[node2]);
//   length = edgeVector[0]*edgeVector[0];
//   for(CFuint iDim=1;iDim <nbDim;iDim++){
//     length += edgeVector[iDim]*edgeVector[iDim];
//     }
//
//   CFreal b = sqrt(length);
//
//   edgeVector = *((*cellNodes)[iNode]) - *((*cellNodes)[node1]);
//   length = edgeVector[0]*edgeVector[0];
//   for(CFuint iDim=1;iDim <nbDim;iDim++){
//     length += edgeVector[iDim]*edgeVector[iDim];
//     }
//
//   CFreal c = sqrt(length);
//   CFreal s = 0.5*(a+b+c);
//   CFreal sinA = 2*(sqrt(s*(s-a)*(s-b)*(s-c)))/(b*c) ;
// //std::cout << sinA << std::endl;
//   return sinA;
// }

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshAdapterSpringAnalogy

  } // namespace Numerics

} // namespace COOLFluiD
