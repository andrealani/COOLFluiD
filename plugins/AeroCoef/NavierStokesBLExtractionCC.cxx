#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Environment/DirPaths.hh"
#include "Framework/PathAppender.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/BadValueException.hh"
#include "FiniteVolume/CellCenterFVM.hh"
#include "FiniteVolume/ComputeDiffusiveFlux.hh"
#include "NavierStokes/NavierStokesVarSet.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "AeroCoef/AeroCoefFVM.hh"
#include "AeroCoef/NavierStokesBLExtractionCC.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Numerics::FiniteVolume;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NavierStokesBLExtractionCC,
		      DataProcessingData,
		      AeroCoefFVMModule>
NavierStokesBLExtractionCCProvider("NavierStokesBLExtractionCC");

//////////////////////////////////////////////////////////////////////////////

void NavierStokesBLExtractionCC::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::string >
    ("OutputFileBL","Name of Output File to write the boundary layer profile.");
  options.addConfigOption< std::vector<CFreal> >
    ("ExtractCoord","Coordinates of the point at which to extract BL");
  options.addConfigOption< bool >
    ("extractBL","Flag if to extract boundary layer profile along the profile");
  options.addConfigOption< CFreal >
    ("BLThickness","Maximum Thickness of the BL for extraction");
  options.addConfigOption< CFreal >
    ("Tolerance","Tolerance for finding the initial coord on a face");

}

//////////////////////////////////////////////////////////////////////////////

NavierStokesBLExtractionCC::NavierStokesBLExtractionCC(const std::string& name) :
  NavierStokesSkinFrictionHeatFluxCC(name)
{
  addConfigOptionsTo(this);

  _extractBLalongProfile = false;
  setParameter("extractBL",&_extractBLalongProfile);

  _outputFileBL = "BLProfile.plt";
  setParameter("OutputFileBL",&_outputFileBL);

  _extractCoord = std::vector<CFreal>();
  setParameter("ExtractCoord",&_extractCoord);

  _BLThickness = 1.;
  setParameter("BLThickness",&_BLThickness);

  _tolerance = 0.000001;
  setParameter("Tolerance",&_tolerance);

}

//////////////////////////////////////////////////////////////////////////////

NavierStokesBLExtractionCC::~NavierStokesBLExtractionCC()
{
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesBLExtractionCC::setup()
{
  CFAUTOTRACE;

  NavierStokesSkinFrictionHeatFluxCC::setup();
  
  //For the moment, we can only extract 2D sections
  cf_assert(_extractCoord.size() == 2);

  _initExtract.resize(_extractCoord.size());
  for(CFuint i=0;i < _extractCoord.size(); i++)
  {
    _initExtract[i] = _extractCoord[i];
  }

}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesBLExtractionCC::executeOnTrs()
{
  CFAUTOTRACE;

  _initCoordFound = false;

  if(_extractBLalongProfile) extractBLalongProfile();
  
  NavierStokesSkinFrictionHeatFluxCC::executeOnTrs();

  // const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  // Execute and save file if needed...
  if(_initCoordFound != true)
    {
      CFLog(NOTICE, "The initial point for the boundary extraction could not be found...\n");
      CFLog(NOTICE, "Check the coordinates entered.(you entered " << _initExtract << ")\n");
    }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesBLExtractionCC::computeExtraValues()
{
  CFAUTOTRACE;
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  cf_assert(dim == DIM_2D);
  
  const vector<Node*>& faceNodes = *m_currFace->getNodes();
  const CFuint nbFaceNodes = faceNodes.size();
  
  RealVector minCoord(dim);
  RealVector maxCoord(dim);
  minCoord = MathTools::MathConsts::CFrealMax();
  maxCoord = -MathTools::MathConsts::CFrealMax();
  
  for (CFuint iNode = 0; iNode < nbFaceNodes; ++iNode) {
    for (CFuint iDim = 0; iDim < dim; ++iDim) {
      minCoord[iDim] = min((*faceNodes[iNode])[iDim], minCoord[iDim]);
      maxCoord[iDim] = max((*faceNodes[iNode])[iDim], maxCoord[iDim]);
    }
  }
  
  if( (_initExtract[XX] >= minCoord[XX] - _tolerance) &&
      (_initExtract[XX] <= maxCoord[XX] + _tolerance) &&
      (_initExtract[YY] >= minCoord[YY] - _tolerance) &&
      (_initExtract[YY] <= maxCoord[YY] + _tolerance))
  {
    extractBoundaryLayerProfile(m_currFace);
    _initCoordFound = true;
  }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesBLExtractionCC::extractBoundaryLayerProfile(GeometricEntity* currFace)
{
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  
  // preparation of the output
  boost::filesystem::path file = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path(_outputFileBL);
  file = PathAppender::getInstance().appendAllInfo(file);

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& fout = fhandle->open(file, ios::app);

  fout << "TITLE = Boundary Layer Profile" << "\n";
  fout << "VARIABLES = x y ";
  std::vector<std::string> varNames = _diffVar->getVarNames();
  for(CFuint iName = 0 ; iName < varNames.size(); iName++)
  {
    fout << varNames[iName] << " ";
  }
  fout << "mu rho uOverUe yPlus uPlus" << endl;

  fout << " ZONE T=\"BL@X= " << _initExtract[0] << "\", ZONETYPE=Ordered DATAPACKING=POINT DT=(SINGLE SINGLE ) " << endl;

  CFuint currFaceID = currFace->getID();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  cf_assert(dim == 2);


  DataHandle< RealVector> nstates = socket_nstates.getDataHandle();
  
  Common::SafePtr<GeometricEntityPool<CellTrsGeoBuilder> >
    geoBuilderCell = m_fvmccData->getCellTrsGeoBuilder();
  
  SafePtr<CellTrsGeoBuilder> geoBuilderCellPtr = geoBuilderCell->getGeoBuilder();
  geoBuilderCellPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);

  CellTrsGeoBuilder::GeoData& geoDataCell = geoBuilderCell->getDataGE();
  geoDataCell.trs = MeshDataStack::getActive()->getTrs("InnerCells");

  ///Get the boundary cell
  //First get the inner state
  State* innerState = currFace->getState(0);
  cf_assert(!innerState->isGhost());

  //Get the cell corresponding to the state
  //Build the GeometricEntity
  geoDataCell.idx = innerState->getLocalID();
  GeometricEntity* nextCell = geoBuilderCell->buildGE();

  //this is a length that should be big enough to span the whole BL but small enough to keep accuracy
  const CFreal vectorLength = 10000000.;
  const CFreal x1 = _initExtract[0];
  const CFreal x2 = _initExtract[0] + m_unitNormal[0] * vectorLength;
  const CFreal y1 = _initExtract[1];
  const CFreal y2 = _initExtract[1] + m_unitNormal[1] * vectorLength;

  // Compute wall state
  std::vector<Node*>* faceNodes = currFace->getNodes();
  cf_assert(faceNodes->size() == 2);

  RealVector vector1 = (*((*faceNodes)[0])) - _initExtract;
  RealVector vector2 = (*((*faceNodes)[1])) - (*((*faceNodes)[0]));
  CFreal dist1 = vector1.norm2();
  CFreal dist2 = vector2.norm2();
  CFreal ratio = dist1/dist2;

  //Compute the values at the intersection
  RealVector extractedValues = nstates[((*faceNodes)[0])->getLocalID()] + ratio*(nstates[((*faceNodes)[1])->getLocalID()] - nstates[((*faceNodes)[0])->getLocalID()]);

  ///@todo this is dangerous
  ///here we assume that diffusive Var and gradient var are the same!

  ///We are at the wall so distance = 0.
  _diffVar->setWallDistance(0.);

  //at the wall the gradient is null
  for(CFuint iEq=0; iEq<nbEqs; iEq++)
  {
    for(CFuint iDim=0; iDim <nbDim; iDim++)
    {
      (*(_gradients[iEq]))[iDim] = 0.;
    }
  }
  const CFreal mu = _diffVar->getDynViscosity(extractedValues, _gradients);
  const CFreal rho = _diffVar->getDensity(extractedValues);

  _muWall = mu;
  _rhoWall = rho;

  fout << _initExtract[0] << " " << _initExtract[1] << " " << extractedValues << " ";
  fout << mu << " ";
  fout << rho << " ";
  fout << " 0. 0. 0." <<"\n";

  //Compute wall static pressure
  ///@todo this is dangerous: here we will assume we are in PUVT
  const CFreal pwall = extractedValues[0];

  // Compute boundary layer edge values (to know when to stop the extraction)
  const CFreal gamma = m_updateVarSet->getModel()->getGamma();
  const CFreal gm1 = gamma - 1.;
  
  const CFreal R = m_updateVarSet->getModel()->getR();
  // unused//  const CFreal rhoInf = m_pInf/(R * _Tinf);
  const CFreal aInf = sqrt(gamma * R * m_TInf);
  const CFreal Minf = m_uInf / aInf;
  const CFreal beta_inf = 1.0 + 0.5 * gm1 * Minf * Minf ;
  const CFreal pt_inf = m_pInf * pow(beta_inf,(gamma/gm1));
  const CFreal Mach_e2 = ((pow(pt_inf/pwall,gm1/gamma)-1.0)*2.0)/gm1;

  // Near the stagnation point, the constant pressure assumption
  // can cause static pressure to be bigger than total pressure!
  CFreal Mach_e = 0.;
  if (Mach_e2 > 0.0)
    Mach_e = sqrt(Mach_e2) ;
  else
    Mach_e = 1.e-10 ;

  const CFreal beta_e = 1.0 + 0.5 * gm1 * Mach_e * Mach_e ;
  const CFreal a_e    = sqrt(beta_inf/beta_e)*aInf;
  const CFreal U_e    = Mach_e * a_e;

  std::cout << "The computed external value for the BL extraction is Ue: " << U_e << std::endl;

  ///loop until BL is over
  bool isEndOfBL(false);
  RealVector intersection = _initExtract;
  RealVector pastIntersection(dim);
  RealVector pastExtractedValues(nbEqs);
  while(!isEndOfBL)
  {
    const vector<State*>* const states = nextCell->getStates();
    
    // all elements in FVM should have only one state
    cf_assert(states->size() == 1);

    const State *const currState = (*states)[0];
    const GeomEntList *const faces = nextCell->getNeighborGeos();
    const CFuint nbFacesInCell = faces->size();

    bool faceFound = false;
    for (CFuint iFace = 0; iFace < nbFacesInCell; ++iFace) {

      const CFuint faceID = (*faces)[iFace]->getID();
      State *const leftState = (*faces)[iFace]->getState(0);
      State *const rightState = (*faces)[iFace]->getState(1);
      std::vector<Node*>* faceNodes = (*faces)[iFace]->getNodes();

      cf_assert(faceNodes->size() == 2);
      const CFreal x3 = (*((*faceNodes)[0]))[0];
      const CFreal x4 = (*((*faceNodes)[1]))[0];
      const CFreal y3 = (*((*faceNodes)[0]))[1];
      const CFreal y4 = (*((*faceNodes)[1]))[1];
//unused//      const CFreal ua = ((x4-x3)*(y1-y3) - (y4-y3)*(x1-x3))/((y4-y3)*(x2-x1) - (x4-x3)*(y2-y1));
      if(((y4-y3)*(x2-x1) - (x4-x3)*(y2-y1)) != 0.){
      const CFreal ub = ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3))/((y4-y3)*(x2-x1) - (x4-x3)*(y2-y1));
// CFout << "ub: " << ub<<"\n";
// CFout << "ubdenom: " << ((y4-y3)*(x2-x1) - (x4-x3)*(y2-y1))<<"\n";
      if((ub <= 1.) && (ub >= 0.) && (faceID != currFaceID))
      {
        faceFound = true;
		pastIntersection = intersection;

        intersection[0] = x3 + ub*(x4-x3);
        intersection[1] = y3 + ub*(y4-y3);

        //temporary stop condition
        RealVector distVector = intersection - _initExtract;
        CFreal dist = distVector.norm2();

        ///Extract the variables at the given point
        //get the distance to the nodal values
        vector1 = (*((*faceNodes)[0])) - intersection;
        vector2 = (*((*faceNodes)[1])) - (*((*faceNodes)[0]));
        dist1 = vector1.norm2();
        dist2 = vector2.norm2();
        ratio = dist1/dist2;

		pastExtractedValues = extractedValues;

        //Compute the values at the intersection
        extractedValues = nstates[((*faceNodes)[0])->getLocalID()] + ratio*(nstates[((*faceNodes)[1])->getLocalID()] - nstates[((*faceNodes)[0])->getLocalID()]);

        //print the values to file
        fout << intersection[0] << " " << intersection[1] << " " ;
        fout << extractedValues;

		//we compute the gradients approximatively
		///@todo improve this
		for(CFuint iEq=0; iEq<nbEqs; iEq++)
		{
           	const CFreal dV = (extractedValues[iEq] - pastExtractedValues[iEq]);
	  		for(CFuint iDim=0; iDim <nbDim; iDim++)
	  		{
	    		const CFreal dX = (intersection[iDim] - pastIntersection[iDim]);
	    		if(dX>MathTools::MathConsts::CFrealEps()) (*(_gradients[iEq]))[iDim] = dV/dX;
	    		else (*(_gradients[iEq]))[iDim] = 0.;
	  		}
		}
        _diffVar->setWallDistance(dist);

        ///@todo this is dangerous
        ///here we assume that diffusive Var and gradient var are the same!
        const CFreal mu = _diffVar->getDynViscosity(extractedValues, _gradients);
        const CFreal rho = _diffVar->getDensity(extractedValues);
        const CFreal U = sqrt(extractedValues[1]*extractedValues[1] + extractedValues[2]*extractedValues[2]);
        fout << mu << " ";
        fout << rho << " ";
        fout << U/U_e << " ";

        //Fixed distance stop condition
        if(dist > _BLThickness) isEndOfBL = true;
        //new BL stop condition
        if(U > 0.99*U_e) isEndOfBL = true;

        CFreal uPlus = sqrt(_tau/_rhoWall);
        fout << " " << sqrt(_rhoWall) * sqrt(_tau) * dist / _muWall << " " << U/uPlus << endl;

        ///Get the next cell ID
        CFuint nextCellID = 0;
        if(leftState == currState)
        {
          nextCellID = rightState->getLocalID();
          if(rightState->isGhost()) isEndOfBL = true;
        }
        else
        {
          nextCellID = leftState->getLocalID();
          if(leftState->isGhost()) isEndOfBL = true;
        }

        ///Build the next cell
        //first release the entity
        geoBuilderCell->releaseGE();
        //then build the new one
        geoDataCell.idx = nextCellID;
        nextCell = geoBuilderCell->buildGE();
        iFace = nbFacesInCell-1;
        currFaceID = faceID;
      }

      }
    }

    cf_assert(faceFound);

    geoBuilderCell->releaseGE();
  }
  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesBLExtractionCC::findStagnationPoint()
{
  CFAUTOTRACE;

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle< RealVector> nstates = socket_nstates.getDataHandle();

  ///Get nodal values on the TRS
  SafePtr<TopologicalRegionSet> currTrs = getCurrentTRS();

  Common::SafePtr<std::vector<CFuint> > trsNodes = currTrs->getNodesInTrs();
  const CFuint nbNodes = trsNodes->size();

  //Loop over the nodal values to find
  //Stagnation point (node with maximum pressure).
  CFreal maxCp = -MathTools::MathConsts::CFrealMax();
  CFreal maxXCoord = -MathTools::MathConsts::CFrealMax();
  RealVector stagnationPointCoord(PhysicalModelStack::getActive()->getDim());
  RealVector trailingEdgeCoord(PhysicalModelStack::getActive()->getDim());

  for (CFuint iNode = 0; iNode < nbNodes; iNode++)
  {
    Node* currNode = nodes[(*trsNodes)[iNode]];
    RealVector currState = nstates[(*trsNodes)[iNode]];

    // unused // const CFreal p = currState[0];
    const CFreal pDim = currState[0] * (m_updateVarSet->getModel()->getPressRef());
    const CFreal rhoInf = m_pInf / (m_updateVarSet->getModel()->getRdim() * m_TInf);
    const CFreal localCp = (pDim - m_pInf)/(0.5 * rhoInf * m_uInf * m_uInf);
    if(localCp > maxCp )
    {
      // No interpolation done, stagnation point is defined at a node.
      maxCp = localCp;
      stagnationPointCoord = *currNode;
      _stagnationPointIdx = iNode;
    }

    if((*currNode)[0] > maxXCoord)
    {
      maxXCoord = (*currNode)[0];
      trailingEdgeCoord = *currNode;
      _trailingEdgeIdx = iNode;
    }

  }

  std::cout << "The stagnation point was found at: " << stagnationPointCoord << std::endl;
  std::cout << "The trailing edge was found at: " << trailingEdgeCoord << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesBLExtractionCC::extractBLalongProfile()
{
  CFAUTOTRACE;

  DataHandle< CFreal> normals = socket_normals.getDataHandle();

  //build the node-face connectivity for the TRS
  buildNodeFaceConnectivity();
  findStagnationPoint();

  //start at the stagnation point
  //find the faces to which the stagnation point belongs
  //if only one, then we have a one face 'object' (flatplate)
  //if two, then we have a two sided 'object' (airfoil)
  const CFuint nbSides = _nodeToFaceConnectivity[_stagnationPointIdx].size();

  Common::SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> >
    geoBuilder = m_fvmccData->getFaceTrsGeoBuilder();

  SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);

  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
  SafePtr<TopologicalRegionSet> currTrs = getCurrentTRS();

  geoData.trs = currTrs;
  geoData.isBFace = true;

  Common::SafePtr<std::vector<CFuint> > trsNodes = currTrs->getNodesInTrs();
  // unused // const CFuint nbNodes = trsNodes->size();
  CFuint currNodeID = (*trsNodes)[_stagnationPointIdx];
  CFuint currNodeIdx = _stagnationPointIdx;
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  //for each side, loop over the TRS from the stagnation point to the trailing edge
  for(CFuint iSide = 0; iSide < nbSides; iSide++)
  {
    CFuint faceIdx = (_nodeToFaceConnectivity[currNodeIdx])[iSide];
    _distanceToStagnation = 0.;
    bool endOfSide(false);
    while(!endOfSide)
    {
      // build the GeometricEntity
      geoData.idx = faceIdx;
      GeometricEntity* const currFace = geoBuilder->buildGE();

      //1. extract the boundary layer profiles perpendicularly to the faces
      std::vector<Node*>* faceNodes = currFace->getNodes();

      _initExtract = 0.5 *(  *((*faceNodes)[0]) +  *((*faceNodes)[1]));

      const CFuint faceID = currFace->getID();
      const CFuint startID = faceID*dim;
      for (CFuint iDim = 0; iDim < dim; ++iDim) {
        m_unitNormal[iDim] = normals[startID + iDim];
      }

      const CFreal faceLength = m_unitNormal.norm2();
      // consider the normal pointing inward the domain
      m_unitNormal *= (-1./faceLength);

      ///Compute the distance to the stagnation point
      _distanceToStagnation += (0.5 * faceLength);

      ///Extract the BL at the center of the face
      extractBoundaryLayerProfile(currFace);
      std::cout << "Distance to stagnation: " << _distanceToStagnation <<std::endl;
      //add the second half of the face to the computation of the distance
      _distanceToStagnation += (0.5 * faceLength);

      //2. Go to the next face
      //Get idx of next node
      for(CFuint iNode=0; iNode < faceNodes->size(); iNode++)
      {
        CFuint nodeID = (*faceNodes)[iNode]->getLocalID();
        if(nodeID != currNodeID)
        {
          currNodeID = nodeID;
          currNodeIdx = _nodesLocalID2TrsIdx.find(nodeID);
          iNode = faceNodes->size();
        }
      }
      if(currNodeIdx == _trailingEdgeIdx){
        endOfSide = true;
      }
      else{
        const CFuint nbFaces = _nodeToFaceConnectivity[currNodeIdx].size();
        for(CFuint iFace = 0; iFace < nbFaces; iFace++)
        {
          //Get idx of next face
          CFuint otherFaceIdx = (_nodeToFaceConnectivity[currNodeIdx])[iFace];
          if(otherFaceIdx != faceIdx)
          {
            faceIdx = otherFaceIdx;
            iFace = nbFaces;
          }
        }
      }

      //release the GeometricEntity
      geoBuilder->releaseGE();
    }
  }

}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesBLExtractionCC::buildNodeFaceConnectivity()
{

  CFAUTOTRACE;

  SafePtr<TopologicalRegionSet> currTrs = getCurrentTRS();
  Common::SafePtr< vector<CFuint> > const trsNodes = currTrs->getNodesInTrs();
  const CFuint nbTrsNodes = trsNodes->size();

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  ///Build the nodes LocalID to TrsIdx map.
  _nodesLocalID2TrsIdx.reserve(nbTrsNodes);
  // fill the map
  for(CFuint iNode = 0; iNode < nbTrsNodes; ++iNode) {
    Node *const currNode = nodes[(*trsNodes)[iNode]];
    const CFuint nodeID = currNode->getLocalID();
    _nodesLocalID2TrsIdx.insert(nodeID,iNode);
  }
  _nodesLocalID2TrsIdx.sortKeys();


  /// resize the socket of pointers to GeomEntity
  _nodeToFaceConnectivity.resize(trsNodes->size());

  Common::SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> >
    geoBuilder = m_fvmccData->getFaceTrsGeoBuilder();

  SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);

  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  geoData.trs = currTrs;
  geoData.isBFace = true;

  CFuint idx = 0;
  for (CFuint iNode = 0; iNode < trsNodes->size(); ++iNode)
  {
    Node *const currNode = nodes[(*trsNodes)[iNode]];

    const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
    for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {

      // build the GeometricEntity
      geoData.idx = iFace;
      GeometricEntity* const currFace = geoBuilder->buildGE();

      bool contained = currFace->containNode(currNode);
      if(contained)
      {
        _nodeToFaceConnectivity[idx].push_back(iFace);
      }
    //release the GeometricEntity
    geoBuilder->releaseGE();
    }

    idx++ ;
  }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoef

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////




