// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Framework/DataProcessing.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "MeshTools/MeshTools.hh"
#include "MeshTools/ComputeWallDistanceNewton.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MeshTools {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<ComputeWallDistanceNewton, DataProcessingData, MeshToolsModule> computeWallDistanceNewtonProvider("ComputeWallDistanceNewton");

//////////////////////////////////////////////////////////////////////////////

void ComputeWallDistanceNewton::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Tolerance","Precision of the distance computation");
  options.addConfigOption< CFreal >("RelaxationFactor","Relaxation Factor for the newton iterations");
  options.addConfigOption< CFuint >("MaxIter","Maximum number of iterations for the newton procedure");
}

//////////////////////////////////////////////////////////////////////////////

ComputeWallDistanceNewton::ComputeWallDistanceNewton(const std::string& name) :
  ComputeWallDistance(name)
{
   addConfigOptionsTo(this);

  _tolerance = 0.0000000001;
   setParameter("Tolerance",&_tolerance);

  _relaxFactor = 1000.;
   setParameter("RelaxationFactor",&_relaxFactor);

  _maxIter = 100;
  setParameter("MaxIter",&_maxIter);

}

//////////////////////////////////////////////////////////////////////////////

ComputeWallDistanceNewton::~ComputeWallDistanceNewton()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeWallDistanceNewton::execute()
{
  CFAUTOTRACE;

  CFLog(INFO, "Computing wall distances using Newton's method...\n");
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle< CFreal> wallDistance = socket_wallDistance.getDataHandle();

  _shortestDistanceComputer.setDimension(PhysicalModelStack::getActive()->getDim()-1);
  _shortestDistanceComputer.setRelaxationFactor(_relaxFactor);
  _shortestDistanceComputer.setTolerance(_tolerance);
  _shortestDistanceComputer.setMaxIter(_maxIter);

  const CFuint nbStates = states.size();
  CFreal minimumDistance;
  CFreal distance;
  for (CFuint iState = 0; iState < nbStates; ++iState)
  {
    minimumDistance = MathTools::MathConsts::CFrealMax();
    for(CFuint iTRS = 0; iTRS < _boundaryTRS.size() ; ++iTRS)
    {
      Common::SafePtr<Framework::TopologicalRegionSet> const faces =
         MeshDataStack::getActive()->getTrs(_boundaryTRS[iTRS]);
      Common::SafePtr<std::vector<CFuint> > nodesInTrs = faces->getNodesInTrs();

      Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
      geoBuilder = getMethodData().getStdTrsGeoBuilder();

      StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
      geoData.trs = faces;
      const CFuint nbFaces = faces->getLocalNbGeoEnts();

      bool found(false);
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
        CFLogDebugMed( "Computing iFace = " << iFace << "\n");

        // build the GeometricEntity
        geoData.idx = iFace;
        GeometricEntity& currFace = *geoBuilder->buildGE();
        vector<Node*>* faceNodes = currFace.getNodes();
        const CFuint nbFaceNodes = faceNodes->size();
        bool faceSelected(false);

        for (CFuint iFaceNode = 0; iFaceNode < nbFaceNodes; ++iFaceNode) {
          _tmpVector = *(nodes[(*faceNodes)[iFaceNode]->getLocalID()]);
          _tmpVector -= states[iState]->getCoordinates();
          distance = _tmpVector.norm2();

          if(distance < 1.1*minimumDistance) faceSelected = true;
          if(distance == 0.){
             minimumDistance = 0.;
             found = true;
             }
          }

        //if the nodes are very close but not the same node
        if((faceSelected) && (found != true)){
          _shortestDistanceComputer.compute(&currFace, states[iState]->getCoordinates(), distance);
          minimumDistance = min(distance, minimumDistance);
        }
        geoBuilder->releaseGE();
      }
    }

    wallDistance[iState] = minimumDistance;
  }

  printToFile();
  CFLog(INFO, "Wall distances computation finished...\n");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace MeshTools

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
