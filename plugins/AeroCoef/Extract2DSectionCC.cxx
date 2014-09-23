#include "Common/PE.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/BadValueException.hh"
#include "Framework/PathAppender.hh"
#include "Environment/DirPaths.hh"
#include "NavierStokes/EulerVarSet.hh"
#include "NavierStokes/NavierStokesVarSet.hh"
#include "FiniteVolume/CellCenterFVM.hh"
#include "FiniteVolume/ComputeDiffusiveFlux.hh"
#include "AeroCoef/AeroCoefFVM.hh"
#include "AeroCoef/Extract2DSectionCC.hh"

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

MethodCommandProvider<Extract2DSectionCC,
		      DataProcessingData,
		      AeroCoefFVMModule>
Extract2DSectionCCProvider("Extract2DSectionCC");

//////////////////////////////////////////////////////////////////////////////

void Extract2DSectionCC::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFuint >
    ("SaveRate","Rate for saving the extracted section.");
  options.addConfigOption< bool >
    ("AppendIter","Append the iteration number to the name of the output file.");
  options.addConfigOption< bool >
    ("AppendTime","Append the time to the name of the output file.");
  options.addConfigOption< std::string >
    ("OutputFile","Name of Output File to write the results.");
  options.addConfigOption< std::vector<CFreal> >
    ("ExtractCoord","Coordinates of the point at which to extract the section");
  options.addConfigOption< bool >
    ("ExtractAlongNormal","Flag if the sectiuon is to be taken normal to the TRS");
  options.addConfigOption< std::vector<CFreal>  >
    ("SectionDirection","If section extracted not normal to the TRS, direction vector of the section");
  options.addConfigOption< CFreal>
    ("Tolerance","Tolerance for searching the initial section point");

}

//////////////////////////////////////////////////////////////////////////////

Extract2DSectionCC::Extract2DSectionCC(const std::string& name) :
  DataProcessingCom(name),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_gstates("gstates"),
  socket_nstates("nstates"),
  socket_normals("normals"),
  _updateVar(CFNULL)
{
  addConfigOptionsTo(this);

  _saveRate = 1;
  setParameter("SaveRate",&_saveRate);

  _appendIter = false;
  setParameter("AppendIter",&_appendIter);

  _appendTime = false;
  setParameter("AppendTime",&_appendTime);

  _outputFile = "extracted2DSection.plt";
  setParameter("OutputFile",&_outputFile);

  _extractCoord = std::vector<CFreal>();
  setParameter("ExtractCoord",&_extractCoord);

  _extractAlongNormal = true;
  setParameter("ExtractAlongNormal",&_extractAlongNormal);

  _sectionDirection = std::vector<CFreal>();
  setParameter("SectionDirection",&_sectionDirection);

  _tolerance = MathTools::MathConsts::CFrealEps();
  setParameter("Tolerance",&_tolerance);

}

//////////////////////////////////////////////////////////////////////////////

Extract2DSectionCC::~Extract2DSectionCC()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
Extract2DSectionCC::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_gstates);
  result.push_back(&socket_nstates);
  result.push_back(&socket_normals);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
Extract2DSectionCC::providesSockets()
{
  // check this
  vector<Common::SafePtr<BaseDataSocketSource> > result;
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void Extract2DSectionCC::setup()
{
  CFAUTOTRACE;

  // suppose that just one space method is available
  SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  SafePtr<Numerics::FiniteVolume::CellCenterFVM> fvmcc = spaceMethod.d_castTo<Numerics::FiniteVolume::CellCenterFVM>();
  cf_assert(fvmcc.isNotNull());

  _fvmccData = fvmcc->getData();

  _updateVar = _fvmccData->getUpdateVar();
  cf_assert(_updateVar.isNotNull());

  _sectionNormal.resize(PhysicalModelStack::getActive()->getDim());

  cf_assert(_extractCoord.size() == 2);

  _initExtract.resize(_extractCoord.size());
  for(CFuint i=0;i < _extractCoord.size(); i++)
  {
    _initExtract[i] = _extractCoord[i];
  }

  if(!_extractAlongNormal){
    cf_assert(_sectionDirection.size() == 2);
    for(CFuint i=0;i < _sectionDirection.size(); i++)
    {
      _sectionNormal[i] = _sectionDirection[i];
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void Extract2DSectionCC::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void Extract2DSectionCC::executeOnTrs()
{
  CFAUTOTRACE;

  CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
  // Execute and save file if needed...
  if (!(iter % _saveRate)) {
    computeValues();
  }
}

//////////////////////////////////////////////////////////////////////////////

void Extract2DSectionCC::computeValues()
{
  CFAUTOTRACE;

  using namespace boost::filesystem;

//   // preparation of the output
//   path fpath = Environment::DirPaths::getInstance().getResultsDir() / _outputFile;
//   fpath = Framework::PathAppender::getInstance().appendParallel( fpath );
//
//   SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
//   ofstream& fout = fhandle->open(fpath);

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle< RealVector> nstates = socket_nstates.getDataHandle();
  DataHandle< CFreal> normals = socket_normals.getDataHandle();

  Common::SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> >
    geoBuilder = _fvmccData->getFaceTrsGeoBuilder();

  SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
  geoBuilderPtr->setDataSockets(socket_states, socket_gstates, socket_nodes);

  FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  SafePtr<TopologicalRegionSet> currTrs = getCurrentTRS();
  CFLog(VERBOSE, "Extract2DSectionCC::computeValues() on TRS "
	<< currTrs->getName() << "\n");

  geoData.trs = currTrs;
  geoData.isBFace = true;

  bool faceFound = false;
  const CFuint nbTrsFaces = currTrs->getLocalNbGeoEnts();
  for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
    CFLogDebugMed("iFace = " << iFace << "\n");

    // build the GeometricEntity
    geoData.idx = iFace;
    GeometricEntity* const currFace = geoBuilder->buildGE();

    const vector<Node*>& faceNodes = *currFace->getNodes();
    const CFuint nbFaceNodes = faceNodes.size();

    CFreal minCoordX = MathTools::MathConsts::CFrealMax();
    CFreal maxCoordX = -MathTools::MathConsts::CFrealMax();
    CFreal minCoordY = MathTools::MathConsts::CFrealMax();
    CFreal maxCoordY = -MathTools::MathConsts::CFrealMax();

    for (CFuint iNode = 0; iNode < nbFaceNodes; ++iNode) {
      minCoordX = min((*faceNodes[iNode])[0], minCoordX);
      minCoordY = min((*faceNodes[iNode])[1], minCoordY);
      maxCoordX = max((*faceNodes[iNode])[0], maxCoordX);
      maxCoordY = max((*faceNodes[iNode])[1], maxCoordY);
    }
/*CFout << "-----------------------------------\n";
CFout << minCoordX << " " <<minCoordY <<"\n";
CFout << maxCoordX << " " <<maxCoordY <<"\n";*/
    if((_initExtract[0] >= minCoordX-MathTools::MathConsts::CFrealEps()) && (_initExtract[0] <= maxCoordX) && (_initExtract[1] >= minCoordY) && (_initExtract[1] <= maxCoordY))
    {
      faceFound = true;
      if(_extractAlongNormal){
        // compute the unit face normal
        const CFuint faceID = currFace->getID();
        const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
        for (CFuint iDim = 0; iDim < PhysicalModelStack::getActive()->getDim(); ++iDim) {
          _sectionNormal[iDim] = normals[startID + iDim];
        }

        // consider the normal pointing inward the domain
        _sectionNormal *= (-1./_sectionNormal.norm2());
      }

      extractSection(currFace);
    }

    geoBuilder->releaseGE();
  }

  if(!faceFound){
    CFLog(NOTICE, "The initial point for the boundary extraction could not be found...\n");
    CFLog(NOTICE, "Check the coordinates. (You entered " << _initExtract << ")\n");
  }

//   fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void Extract2DSectionCC::extractSection(GeometricEntity* currFace)
{
  CFAUTOTRACE;

  using namespace boost::filesystem;

  // preparation of the output
  path fpath = Environment::DirPaths::getInstance().getResultsDir() / path( _outputFile );
  fpath = Framework::PathAppender::getInstance().appendAllInfo( fpath , _appendIter, _appendTime);

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& fout = fhandle->open(fpath,ios::app);

  fout << "TITLE = SectionExtracted" << "\n";
  fout << "VARIABLES = x y ";

  std::vector<std::string> varNames = _updateVar->getVarNames();
  for(CFuint iName = 0 ; iName < varNames.size(); iName++)
  {
    fout << varNames[iName] << " ";
  }
  fout << "\n";

  CFuint currFaceID = currFace->getID();
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  cf_assert(dim == 2);

  DataHandle< RealVector> nstates = socket_nstates.getDataHandle();

  Common::SafePtr<GeometricEntityPool<CellTrsGeoBuilder> >
    geoBuilderCell = _fvmccData->getCellTrsGeoBuilder();

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
  const CFreal vectorLength = 1000000.;
  const CFreal x1 = _initExtract[0];
  const CFreal x2 = _initExtract[0] + _sectionNormal[0] * vectorLength;
  const CFreal y1 = _initExtract[1];
  const CFreal y2 = _initExtract[1] + _sectionNormal[1] * vectorLength;

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

  fout << _initExtract[0] << " " << _initExtract[1] << " " << extractedValues <<"\n";

  ///loop until mesh is over
  bool isEndOfSection(false);
  RealVector intersection(dim);
  while(!isEndOfSection)
  {
    const vector<State*>* const states = nextCell->getStates();
//unused//    const CFuint cellID = nextCell->getID();

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
      const CFreal ub = ((x2-x1)*(y1-y3) - (y2-y1)*(x1-x3))/((y4-y3)*(x2-x1) - (x4-x3)*(y2-y1));

      if((ub <= 1.) && (ub >= 0.) && (faceID != currFaceID))
      {
        faceFound = true;

        intersection[0] = x3 + ub*(x4-x3);
        intersection[1] = y3 + ub*(y4-y3);

        ///Extract the variables at the given point
        //get the distance to the nodal values
        vector1 = (*((*faceNodes)[0])) - intersection;
        vector2 = (*((*faceNodes)[1])) - (*((*faceNodes)[0]));
        dist1 = vector1.norm2();
        dist2 = vector2.norm2();
        ratio = dist1/dist2;

        //Compute the values at the intersection
        extractedValues = nstates[((*faceNodes)[0])->getLocalID()] + ratio*(nstates[((*faceNodes)[1])->getLocalID()] - nstates[((*faceNodes)[0])->getLocalID()]);

        //print the values to file
        fout << intersection[0] << " " << intersection[1] << " " << extractedValues << "\n";

        ///Get the next cell ID
        CFuint nextCellID = 0;
        if(leftState == currState)
        {
          if(!rightState->isGhost()) nextCellID = rightState->getLocalID();
          else isEndOfSection = true;
        }
        else
        {
          if(!leftState->isGhost()) nextCellID = leftState->getLocalID();
          else isEndOfSection = true;
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

    cf_assert(faceFound);

    geoBuilderCell->releaseGE();
  }
  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void Extract2DSectionCC::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////




