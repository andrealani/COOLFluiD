#include <boost/progress.hpp>

#include "Common/FilesystemException.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/DataHandle.hh"
#include "Framework/SubSystemStatus.hh"
#include "Environment/DirPaths.hh"
#include "Framework/MeshData.hh"
#include "Framework/CommandGroup.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "SubSystemCoupler/StdMeshMatcherWrite2.hh"
#include "SubSystemCoupler/SubSystemCoupler.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdMeshMatcherWrite2, SubSysCouplerData, SubSystemCouplerModule> StdMeshMatcherWrite2Provider("StdMeshMatcherWrite2");

//////////////////////////////////////////////////////////////////////////////

void StdMeshMatcherWrite2::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< bool> ("MatchModifiedMesh","Match mesh moved by the value of the states.");
}

//////////////////////////////////////////////////////////////////////////////

StdMeshMatcherWrite2::StdMeshMatcherWrite2(const std::string& name) :
  CouplerCom(name),
  _sockets(),
  _matchingFace(static_cast<Framework::TopologicalRegionSet*>(CFNULL),CFNULL),
  _shapeFunctionAtCoord()
{
   addConfigOptionsTo(this);

  _shiftStates = false;
   setParameter("MatchModifiedMesh",&_shiftStates);
}

//////////////////////////////////////////////////////////////////////////////

StdMeshMatcherWrite2::~StdMeshMatcherWrite2()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdMeshMatcherWrite2::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result = _sockets.getAllSinkSockets();;

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdMeshMatcherWrite2::configure ( Config::ConfigArgs& args )
{
  CouplerCom::configure(args);

  typedef std::vector<Common::SafePtr<CommandGroup> > InterfaceList;
  InterfaceList interfaces = getMethodData().getInterfaces();

  try {
    const std::string nsp = getMethodData().getNamespace();
  InterfaceList::iterator itr = interfaces.begin();
  for(; itr != interfaces.end(); ++itr) {

    const std::vector<std::string>& comNames = (*itr)->getComNames();

    // check if this command applies to this Interface
    if(count(comNames.begin(),comNames.end(),getName()) != 0) {

      const std::string interfaceName = (*itr)->getName();
      
      for (CFuint iProc = 0; iProc < Common::PE::GetPE().GetProcessorCount(nsp); ++iProc)
      {
        const vector<std::string> otherTrsNames = getMethodData().getCoupledSubSystemsTRSNames(interfaceName);
        for (CFuint iTRS=0; iTRS< otherTrsNames.size(); ++iTRS)
        {
          /// Create the datahandle for the otherSubSystem Coordinates
          const std::string otherTrsName = otherTrsNames[iTRS];

          /// Get the datahandle for the coordinates of the OTHER SubSystem
          vector<std::string> socketCoordNames = getMethodData().getOtherCoupledCoordName(interfaceName,otherTrsName,iProc);
          vector<std::string> socketAcceptNames = getMethodData().getOtherCoupledAcceptedName(interfaceName,otherTrsName,iProc);
          vector<std::string> socketDataNames = getMethodData().getOtherCoupledDataName(interfaceName,otherTrsName,iProc);
          for(CFuint iType=0;iType < socketCoordNames.size();iType++)
          {
            _sockets.createSocketSink<RealVector>(socketCoordNames[iType]);
            _sockets.createSocketSink<CFreal>(socketAcceptNames[iType]);
            _sockets.createSocketSink<RealVector>(socketDataNames[iType]);
          }
        } // other trs
      } //loop over processors
    } // if check
  } // interfaces

  } catch (Exception& e)
  {
    CFout << e.what() << "\n" << CFendl;
    throw;
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdMeshMatcherWrite2::execute()
{
  CFAUTOTRACE;
  
  const std::string nsp = getMethodData().getNamespace();
  for (CFuint iProc = 0; iProc < Common::PE::GetPE().GetProcessorCount(nsp); ++iProc) {
    executeWrite(iProc);
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdMeshMatcherWrite2::executeWrite(const CFuint iProc)
{
  CFAUTOTRACE;

  if(_shiftStates) shiftStatesCoord();

  /// Get the names of the interfaces, subsystems
  const std::string interfaceName = getCommandGroupName();
  vector<std::string> otherTrsNames = getMethodData().getCoupledSubSystemsTRSNames(interfaceName);

  /// If the geometry is non matching, then get the threshold for acceptance
  /// else set the threshold to a very high value
  const bool isNonMatchingGeometry = getMethodData().isInterfaceGeometryNonMatching(interfaceName);
  CFreal nonMatchingGeometryThreshold = MathTools::MathConsts::CFrealMax();
  if(isNonMatchingGeometry) nonMatchingGeometryThreshold = getMethodData().getNonMatchingGeometryThreshold(interfaceName);

  /// Get the geometric entity builder
  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  /// Loop over the TRSs of the Interface
  for (CFuint iTRS=0; iTRS < otherTrsNames.size(); iTRS++)
  {
    CFLogInfo ("Matching for TRS [" << otherTrsNames[iTRS] << "]\n");

    const vector<std::string> socketCoordNames =
      getMethodData().getOtherCoupledCoordName(interfaceName,otherTrsNames[iTRS],iProc);

    const vector<std::string> socketAcceptNames =
      getMethodData().getOtherCoupledAcceptedName(interfaceName,otherTrsNames[iTRS],iProc);
    const vector<std::string> socketDataNames =
      getMethodData().getOtherCoupledDataName(interfaceName,otherTrsNames[iTRS],iProc);

    for(CFuint iType=0; iType < socketCoordNames.size();iType++)
    {

      CFLogInfo("Matching [" << socketCoordNames[iType] << "]\n");

      DataHandle< RealVector> interfaceCoords =
        _sockets.getSocketSink<RealVector>(socketCoordNames[iType])->getDataHandle();
      DataHandle< CFreal> isAccepted =
        _sockets.getSocketSink<CFreal>(socketAcceptNames[iType])->getDataHandle();
      DataHandle< RealVector> interfaceData =
        _sockets.getSocketSink<RealVector>(socketDataNames[iType])->getDataHandle();

      // Get the coupledGeoEntities and resize
      const CFuint otherNbStates = interfaceCoords.size();
      SubSysCouplerData::CoupledGeoEntities* coupledGeoEntities =
        getMethodData().getCoupledInterfaces(interfaceName,iTRS, iType, iProc);

      (*coupledGeoEntities).resize(otherNbStates);

      // Counter for the number of rejected states
      CFuint rejectedStates = 0;

      // For each state, find the interpolated value
      // on the other boundary
      boost::progress_display progress (otherNbStates);
      for (CFuint iState = 0; iState < otherNbStates; ++iState, ++progress) {

        CFint stateID = -1;
        CFuint idx = iState-rejectedStates;
        bool accepted = false;
        _minimumDistanceOnFace = MathTools::MathConsts::CFrealMax();
        _minimumDistanceOffFace = MathTools::MathConsts::CFrealMax();

        const RealVector coord = interfaceCoords[iState];
        RealVector coordProj(coord.size());

        // Pair each fluid point with the closest wet structural element
        nodeToElementPairing(coord,stateID,coordProj);

        // if the projection of the point is outside all faces, take the closest node
        if (stateID != -1){
          if(nonMatchingGeometryThreshold > _minimumDistanceOffFace) accepted = true;

          if(accepted)
          {
            isAccepted[iState] = _minimumDistanceOffFace;
            // Resize
            (*coupledGeoEntities)[idx].third.resize(_shapeFunctionAtCoord.size());
            (*coupledGeoEntities)[idx].fourth.resize(PhysicalModelStack::getActive()->getDim());
            // set values
            (*coupledGeoEntities)[idx].first = _matchingFace.first;
            (*coupledGeoEntities)[idx].second = _matchingFace.second;
            (*coupledGeoEntities)[idx].third = _shapeFunctionAtCoord;
            (*coupledGeoEntities)[idx].fourth = coordProj;
            CFout << "Node: "<< coord << " does not project on a face and was projected on coord: " << coordProj << "\n";
            CFLogDebugMin("...Distance: " << _minimumDistanceOffFace << "\n");
          }
          else
          {
            isAccepted[iState] = -1.;
            rejectedStates++;
            CFout << "Node doesn't project on a face and is rejected because distance is: "<< _minimumDistanceOffFace << "\n";
          }
        }
        // else continue
        else {
          if(nonMatchingGeometryThreshold > _minimumDistanceOnFace) accepted = true;

          if(accepted)
          {
            isAccepted[iState] = _minimumDistanceOnFace;
            // Determine the natural coordinates of the projection of the fluid point
            // and compute the shape function at this location
            // Resize
            (*coupledGeoEntities)[idx].third.resize(_shapeFunctionAtCoord.size());
            (*coupledGeoEntities)[idx].fourth.resize(PhysicalModelStack::getActive()->getDim());
            // set values
            (*coupledGeoEntities)[idx].first = _matchingFace.first;
            (*coupledGeoEntities)[idx].second = _matchingFace.second;
            (*coupledGeoEntities)[idx].third = _shapeFunctionAtCoord;
            (*coupledGeoEntities)[idx].fourth = coordProj;

            CFLogDebugMin("Node: "<< coord << " projects on a face at coord: " << coordProj << "\n");
            CFLogDebugMin("Accepted: " << coord << "\n");
            CFLogDebugMin("...and associated to: " << coordProj << "\n");
            CFLogDebugMin("...idx: " << idx << "\n");
            CFLogDebugMin("...Distance: " << _minimumDistanceOnFace << "\n");
            CFLogDebugMin("...ShapeFunction: " << (*coupledGeoEntities)[idx].third << "\n");

            geoData.trs = _matchingFace.first;
            geoData.idx = _matchingFace.second;
            GeometricEntity& matchingFaceGeo = *geoBuilder->buildGE();
            std::vector<State*>* states = matchingFaceGeo.getStates();
            CFLogDebugMin("...nbStates: " << states->size() << "\n");
            CFLogDebugMin("...State[0]: " << (*states)[0]->getCoordinates() << "\n");
            CFLogDebugMin("...State[1]: " << (*states)[1]->getCoordinates() << "\n");
            if(states->size() > 2) { CFLogDebugMin("...State[2]: " << (*states)[2]->getCoordinates() << "\n"); }
            if(states->size() > 3) { CFLogDebugMin("...State[3]: " << (*states)[3]->getCoordinates() << "\n"); }
            if(states->size() > 4) { CFLogDebugMin("...State[4]: " << (*states)[4]->getCoordinates() << "\n"); }
            if(states->size() > 5) { CFLogDebugMin("...State[5]: " << (*states)[5]->getCoordinates() << "\n"); }

//if((coord[0] > 1.17)&&(coord[1] < 0.75)){
if(false){
            CFout <<"Node: "<< coord << " projects on a face at coord: " << coordProj << "\n";
            CFout <<"Accepted: " << coord << "\n";
            CFout <<"...and associated to: " << coordProj << "\n";
            CFout <<"...idx: " << idx << "\n";
            CFout <<"...Distance: " << _minimumDistanceOnFace << "\n";
            CFout <<"...ShapeFunction: " << (*coupledGeoEntities)[idx].third << "\n";
            CFout <<"...nbStates: " << states->size() << "\n";
            CFout <<"...State[0]: " << (*states)[0]->getCoordinates() << "\n";
            CFout <<"...State[1]: " << (*states)[1]->getCoordinates() << "\n";
            if(states->size() > 2) CFout <<"...State[2]: " << (*states)[2]->getCoordinates() << "\n";
            if(states->size() > 3) CFout <<"...State[3]: " << (*states)[3]->getCoordinates() << "\n";
            if(states->size() > 4) CFout <<"...State[4]: " << (*states)[4]->getCoordinates() << "\n";
            if(states->size() > 5) CFout <<"...State[5]: " << (*states)[5]->getCoordinates() << "\n";

            CFout <<"...Node[0]: " << *((*matchingFaceGeo.getNodes())[0]) << "\n";
            CFout <<"...Node[1]: " << *((*matchingFaceGeo.getNodes())[1]) << "\n";
            if(states->size() > 2) CFout <<"...Node[2]: " << *((*matchingFaceGeo.getNodes())[2]) << "\n";

          RealVector mappedCoord = matchingFaceGeo.computeMappedCoordFromCoord(coordProj);
CFout << "mappedCoord: "<< mappedCoord  << "\n";
}

            geoBuilder->releaseGE();
          }
          else
          {
            isAccepted[iState] = -1.;
            rejectedStates++;
            CFLogDebugMin("Node: "<< coord << " projects on a face at coord: " << coordProj << "\n");
            CFout << "Node projects on a face but is rejected because distance is: " << _minimumDistanceOnFace << "\n";
          }
        } //end else
      } //end loop over otherStates

      ///We know the number of accepted states -> can now resize the datahandle
      interfaceData.resize(otherNbStates - rejectedStates);

      ///Write the file with the info isAccepted
      writeIsAcceptedFile(socketAcceptNames[iType]);
    } //end of loop over data transfer coord type
  } // end loop over the OtherTRS

  if(_shiftStates) unshiftStatesCoord();
}

//////////////////////////////////////////////////////////////////////////////

void StdMeshMatcherWrite2::shiftStatesCoord()
{
  /// Loop over all the nodes of all the TRS's of this command
  vector< SafePtr<TopologicalRegionSet> > trs = getTrsList();
  vector< Common::SafePtr<TopologicalRegionSet> >::iterator iTRS;

  ///@todo remove this this is very dirty!!!!!!!
  DataHandle < Framework::State*, Framework::GLOBAL > states =  MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle();
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  for (iTRS = trs.begin(); iTRS != trs.end(); ++iTRS)
  {
    const CFuint nbStates = (*iTRS)->getNbStatesInTrs();
    Common::SafePtr<std::vector<CFuint> > trsStates = (*iTRS)->getStatesInTrs();
    for(CFuint iState = 0; iState < nbStates; ++iState)
    {
      Node* node = states[(*trsStates)[iState]]->getNodePtr();
      for(CFuint iDim = 0; iDim < nbDim; ++iDim)
      {
        (*node)[iDim] = (*node)[iDim] + (*(states[(*trsStates)[iState]]))[iDim];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdMeshMatcherWrite2::unshiftStatesCoord()
{
  /// Loop over all the nodes of all the TRS's of this command
  vector< SafePtr<TopologicalRegionSet> > trs = getTrsList();
  vector< Common::SafePtr<TopologicalRegionSet> >::iterator iTRS;

  ///@todo remove this this is very dirty!!!!!!!
  DataHandle < Framework::State*, Framework::GLOBAL > states =  MeshDataStack::getActive()->getStateDataSocketSink().getDataHandle();

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();

  for (iTRS = trs.begin(); iTRS != trs.end(); ++iTRS)
  {
    const CFuint nbStates = (*iTRS)->getNbStatesInTrs();
    Common::SafePtr<std::vector<CFuint> > trsStates = (*iTRS)->getStatesInTrs();
    for(CFuint iState = 0; iState < nbStates; ++iState)
    {
      Node* node = states[(*trsStates)[iState]]->getNodePtr();
      for(CFuint iDim = 0; iDim < nbDim; ++iDim)
      {
        (*node)[iDim] = (*node)[iDim] -  (*(states[(*trsStates)[iState]]))[iDim];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdMeshMatcherWrite2::nodeToElementPairing(const RealVector& coord, CFint& stateID, RealVector& coord_Proj)
{
  CFAUTOTRACE;

  const CFuint dim = coord.size();
  RealVector normal(dim);
  RealVector x0(dim);
  RealVector x1(dim);
  RealVector v(dim);
  RealVector v1(dim);
  RealVector w(dim);
  RealVector project(dim);
  RealVector project2(dim);
  RealVector projectCoord(dim);

RealVector faceVector(dim);
RealVector nodeStateVector(dim);
RealVector Edge0(dim);
RealVector Edge1(dim);
RealVector vectorD(dim);

//   CFreal tempV, tempW;
  bool isOnFace(false);
  bool projectedOnAFace(false);
  bool found(false);

  /// Loop over all the faces of all the TRS's of this command
  vector< SafePtr<TopologicalRegionSet> > trs = getTrsList();
  vector< SafePtr<TopologicalRegionSet> >::iterator iTRS;

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
CFLogDebugMed("Trying to match coord: " << coord << "\n");
  for (iTRS = trs.begin(); iTRS != trs.end(); ++iTRS)
  {
    const CFuint nbGeos = (*iTRS)->getLocalNbGeoEnts();
    const std::string currentTRSName = (*iTRS)->getName();
    geoData.trs = (*iTRS);
CFLogDebugMed( "Looping on face of TRS: " << currentTRSName << "\n");
    for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt)
    {
      // build the GeometricEntity
      geoData.idx = iGeoEnt;
      GeometricEntity& currFace = *geoBuilder->buildGE();

if(dim == 2){
      /// Check if node defined by coord has his projection in face...
      // 1 - compute the normal to the face
      // compute the face normal in 2D

      x0 = *(currFace.getNode(0));
      x1 = *(currFace.getNode(1));
      v1 = x1-x0;

      normal = currFace.computeAvgCellNormal();
      normal.normalize();

      // 2 - project vector xP-x0 (=v) on the normal from x0 -> w = (n*v)*v
      v = coord - x0;
      w = normal * (MathTools::MathFunctions::innerProd(v,normal));

      // 3 - Obtain the vector xProj-x0 (= v - w)
      project = v-w ;
      projectCoord = project + x0 ;

      isOnFace = currFace.isInElement(projectCoord);

      // 4 - Now we have the point projected on the plane of the face
      //     We have to check if the point is inside the face

      if (isOnFace)
      {
        // Compute the distance of projected point to the nodes and take minimum
        project2 = (projectCoord - coord);
        const CFreal distanceToFace = project2.norm2();

        if ((distanceToFace < _minimumDistanceOnFace) &&(distanceToFace < _minimumDistanceOffFace) ){
          coord_Proj = projectCoord;
          _matchingFace.first = (*iTRS);
          _matchingFace.second = iGeoEnt;

          RealVector mappedCoord = currFace.computeMappedCoordFromCoord(coord_Proj);
          _shapeFunctionAtCoord.resize(currFace.nbStates());
          _shapeFunctionAtCoord = currFace.computeShapeFunctionAtMappedCoord(mappedCoord);
          RealVector geoshapeFunctionAtCoord = currFace.computeGeoShapeFunctionAtMappedCoord(mappedCoord);

/*RealVector coord2 = currFace.computeCoordFromMappedCoord(mappedCoord) - coord_Proj;
CFout << "mappedCoord: "<< mappedCoord  << "\n";
CFout << "_geoshapeFunctionAtCoord: "<< geoshapeFunctionAtCoord  << "\n";
CFout << "_shapeFunctionAtCoord: "<< _shapeFunctionAtCoord  << "\n";
CFout << "distanceToFace: "<< distanceToFace  << "\n";*/
//          CFLogDebugMed("Projects On A Face with Distance: "<< distanceToFace << "\n");
  /*CFout << "Face iGeoEnt: " << iGeoEnt << "\n";
  CFout << "Node 0: " << x0 << "\n";
  CFout << "Node 1: " << x1 << "\n";
  CFout << "Coord Proj: " << coord_Proj << "\n";
  CFout << "_shapeFunctionAtCoord: " << _shapeFunctionAtCoord << "\n";*/
          _minimumDistanceOnFace = distanceToFace;
          //once is has been projected on a face, no need of this
          _minimumDistanceOffFace = -1.;
          projectedOnAFace = true;
          found = true;
          CFLogDebugMed("Has been Projected On A Face - MinDistance: "<< _minimumDistanceOnFace<< "\n");
        }
      }
      else
      {
//         if(!projectedOnAFace){
          // 6 - the point projection might fall outside all faces
          // Compute the coordinates of projected point
          bool compute(false);
          for(CFuint iState=0; iState < currFace.nbStates(); ++iState)
          {
            v = coord - ((currFace.getState(iState)->getCoordinates()));
            const CFreal distance2Node = v.norm2();
/*            if(distance2Node < _minimumDistanceOffFace){

              _matchingFace.first = (*iTRS);
              _matchingFace.second = iGeoEnt;
              _shapeFunctionAtCoord.resize(currFace.nbStates());
              _shapeFunctionAtCoord = 0.;
              _shapeFunctionAtCoord[iState] = 1.;
              stateID = iState;
              _minimumDistanceOffFace = distance2Node;
              coord_Proj = (currFace.getState(iState)->getCoordinates());
              found = true;
              CFLogDebugMed("Not Projected - MinDistance: "<< _minimumDistanceOffFace<< "\n");
            }*/
            if((distance2Node < _minimumDistanceOffFace) && (distance2Node < _minimumDistanceOnFace)) {
              compute=true;
              stateID = iState;
              _minimumDistanceOffFace = distance2Node;
              coord_Proj = (currFace.getState(iState)->getCoordinates());
            }
          }
          if(compute){
            _matchingFace.first = (*iTRS);
            _matchingFace.second = iGeoEnt;
            _shapeFunctionAtCoord.resize(currFace.nbStates());
            _shapeFunctionAtCoord = 0.;
            _shapeFunctionAtCoord[stateID] = 1.;
            found = true;

//             CFreal sum = 0;
//             for(CFuint iState=0; iState < currFace.nbStates(); ++iState)
//             {
//               v = coord - ((currFace.getState(iState)->getCoordinates()));
//               const CFreal distance2Node = v.norm2();
//               _shapeFunctionAtCoord[iState] = 1./distance2Node;
//               sum += 1./distance2Node;
//             }
//             _shapeFunctionAtCoord /= sum;

//           }
        }
      }
}
        if(dim == 3)
        {
          cf_assert(currFace.getNodes()->size() == 3);

          //Taken from: Distance Between Point and Triangle in 3D
          //David Eberly, Geometric Tools, LLC
          //http://www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf

          const Node& firstNode = *currFace.getNode(0);
          Edge0 = *currFace.getNode(1) - *currFace.getNode(0);
          Edge1 = *currFace.getNode(2) - *currFace.getNode(0);
          vectorD = *currFace.getNode(0) - coord;
          CFreal a = MathFunctions::innerProd(Edge0, Edge0);
          CFreal b = MathFunctions::innerProd(Edge0, Edge1);
          CFreal c = MathFunctions::innerProd(Edge1, Edge1);
          CFreal d = MathFunctions::innerProd(Edge0, vectorD);
          CFreal e = MathFunctions::innerProd(Edge1, vectorD);
//          CFreal f = MathFunctions::innerProd(vectorD, vectorD);

          CFreal det = std::abs(a*c-b*b);
          CFreal s = b*e-c*d;
          CFreal t = b*d-a*e;
          CFint region = -1;

          if ( s+t <= det )
          {
            if ( s < 0 ){ if ( t < 0 ) { region = 4; } else { region = 3; } }
            else if ( t < 0 ) { region = 5; }
            else { region = 0; }
          }
          else
          {
            if ( s < 0 ) { region = 2; }
            else if ( t < 0 ) { region = 6; }
            else { region = 1; }
          }

          //Smallest distance point to the face is inside triangle
	  if(region == 0){
            CFreal invDet = 1/det;
            s *= invDet;
            t *= invDet;
          }
	  if(region == 1){
            // F(s) = Q(s,1-s) = (a-2b+c)s^2 + 2(b-c+d-e)s + (c+2e+f)
            // F'(s)/2 = (a-2b+c)s + (b-c+d-e)
            // F'(S) = 0 when S = (c+e-b-d)/(a-2b+c)
            // a-2b+c = |E0-E1|^2 > 0, so only sign of c+e-b-d need be considered
	    CFreal numer = c+e-b-d;

	    if ( numer <= 0 )
	    {
              s = 0;
            }
            else
            {
              CFreal denom = a-2*b+c; // positive quantity
              s = ( numer >= denom ? 1 : numer/denom );
            }
            t = 1-s;
          }
	  if((region == 3) || (region == 5)){
	     // F(t) = Q(0,t) = ct^2 + 2et + f
	     // F'(t)/2 = ct+e
	     // F'(T) = 0 when T = -e/c
	     s = 0;
	     t = ( e >= 0 ? 0 : ( -e >= c ? 1 : -e/c ) );
          }
	  if((region == 2) || (region == 4)|| (region == 6)){
	    // Grad(Q) = 2(as+bt+d,bs+ct+e)
	    // (0,-1)*Grad(Q(0,1)) = (0,-1)*(b+d,c+e) = -(c+e)
	    // (1,-1)*Grad(Q(0,1)) = (1,-1)*(b+d,c+e) = (b+d)-(c+e)
	    // min on edge s+t=1 if (1,-1)*Grad(Q(0,1)) < 0 )
	    // min on edge s=0 otherwise
	    CFreal tmp0 = b+d;
	    CFreal tmp1 = c+e;
	    if ( tmp1 > tmp0 ) // minimum on edge s+t=1
	    {
	      CFreal numer = tmp1 - tmp0;
	      CFreal denom = a-2*b+c;
	      s = ( numer >= denom ? 1 : numer/denom );
	      t = 1-s;
	    }
	    else // minimum on edge s=0
	    {
	      s = 0;
	      t = ( tmp1 <= 0 ? 1 : ( e >= 0 ? 0 : -e/c ) );
  	    }
          }
          cf_assert(region > -1);
//CFout << "Region: " << region << "\n";
          projectCoord = coord - (firstNode + s*Edge0 + t*Edge1);
//CFout << "projectCoord: " << projectCoord << "\n";
	  CFreal stateFaceDistance = projectCoord.norm2();

        if (stateFaceDistance < _minimumDistanceOnFace){
          coord_Proj = (firstNode + s*Edge0 + t*Edge1);
          _matchingFace.first = (*iTRS);
          _matchingFace.second = iGeoEnt;

          RealVector mappedCoord = currFace.computeMappedCoordFromCoord(coord_Proj);
//           mappedCoord[0] = s;
//           mappedCoord[1] = t;
//CFout << "s: " << s << " t: "<< t  << "\n";
          _shapeFunctionAtCoord.resize(currFace.nbStates());
          _shapeFunctionAtCoord = currFace.computeShapeFunctionAtMappedCoord(mappedCoord);
          RealVector geoshapeFunctionAtCoord = currFace.computeGeoShapeFunctionAtMappedCoord(mappedCoord);

/*RealVector coord2 = currFace.computeCoordFromMappedCoord(mappedCoord) - coord_Proj;
CFout << "mappedCoord: "<< mappedCoord  << "\n";
CFout << "_geoshapeFunctionAtCoord: "<< geoshapeFunctionAtCoord  << "\n";
CFout << "_shapeFunctionAtCoord: "<< _shapeFunctionAtCoord  << "\n";
CFout << "distanceToFace: "<< distanceToFace  << "\n";*/
//          CFLogDebugMed("Projects On A Face with Distance: "<< distanceToFace << "\n");
  /*CFout << "Face iGeoEnt: " << iGeoEnt << "\n";
  CFout << "Node 0: " << x0 << "\n";
  CFout << "Node 1: " << x1 << "\n";
  CFout << "Coord Proj: " << coord_Proj << "\n";
  CFout << "_shapeFunctionAtCoord: " << _shapeFunctionAtCoord << "\n";*/
          _minimumDistanceOnFace = stateFaceDistance;
          //once is has been projected on a face, no need of this
          _minimumDistanceOffFace = -1.;
          projectedOnAFace = true;
          found = true;
        }
      }



      //release the GeometricEntity
      geoBuilder->releaseGE();
    } // end of loop over faces
  } // end of loop over TRS's

  if(projectedOnAFace == true){
    stateID = -1;
  }
}

//////////////////////////////////////////////////////////////////////////////

void StdMeshMatcherWrite2::writeIsAcceptedFile(const std::string dataHandleName)
{
  CFAUTOTRACE;

// CFout << "Writing the isAccepted file: " << dataHandleName << "\n";

  DataHandle< CFreal > isAccepted = _sockets.getSocketSink<CFreal>(dataHandleName)->getDataHandle();

  boost::filesystem::path dpath (dataHandleName);
  boost::filesystem::path nameOutputFile =
    Environment::DirPaths::getInstance().getResultsDir() / dpath;

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& fout = fhandle->open(nameOutputFile);

  const CFuint nbStates = isAccepted.size();

  fout << nbStates << "\n";
  for (CFuint iState = 0; iState < nbStates; ++iState) {
    fout << isAccepted[iState] << "\n";
  }

  fhandle->close();
// CFout << "Done Writing the isAccepted file: " << dataHandleName << "\n";
}


//////////////////////////////////////////////////////////////////////////////


    } // namespace MeshMatcher

  } // namespace Numerics

} // namespace COOLFluiD
