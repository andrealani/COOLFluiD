#ifdef CF_HAVE_BOOST_1_85
#define BOOST_TIMER_ENABLE_DEPRECATED
#endif
#include <boost/progress.hpp>

#include "Common/FilesystemException.hh"
#include "Environment/SingleBehaviorFactory.hh"
#include "Environment/FileHandlerInput.hh"
#include "Environment/FileHandlerOutput.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/DataHandle.hh"
#include "Framework/SubSystemStatus.hh"
#include "Environment/DirPaths.hh"
#include "Framework/CommandGroup.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "SubSystemCoupler/SubSystemCoupler.hh"
#include "SubSystemCoupler/NewtonMeshMatcherWrite.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace SubSystemCoupler {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NewtonMeshMatcherWrite, SubSysCouplerData, SubSystemCouplerModule> NewtonMeshMatcherWriteProvider("NewtonMeshMatcherWrite");

//////////////////////////////////////////////////////////////////////////////

void NewtonMeshMatcherWrite::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("MaxAccuracy","Maximum Accuracy for finding the projected point.");
   options.addConfigOption< CFuint >("NbSelectedFaces","Number of preselected faces.");
   options.addConfigOption< CFuint >("MaxNewtonSteps","Maximum Number of Newton Steps for finding the projected point.");
   options.addConfigOption< bool> ("MatchModifiedMesh","Match mesh moved by the value of the states.");
}

//////////////////////////////////////////////////////////////////////////////

NewtonMeshMatcherWrite::NewtonMeshMatcherWrite(const std::string& name) :
  StdMeshMatcherWrite(name),
  _numericalJacob(new NumericalJacobian("NumericalJacobian"))
{
   addConfigOptionsTo(this);

  _maxAccuracy = 1.e-10;
   setParameter("MaxAccuracy",&_maxAccuracy);

  _maxNewtonSteps = 50;
   setParameter("MaxNewtonSteps",&_maxNewtonSteps);

  _nbSelectedFaces = 5;
   setParameter("NbSelectedFaces",&_nbSelectedFaces);

  _shiftStates = false;
   setParameter("MatchModifiedMesh",&_shiftStates);

}

//////////////////////////////////////////////////////////////////////////////

NewtonMeshMatcherWrite::~NewtonMeshMatcherWrite()
{
}

//////////////////////////////////////////////////////////////////////////////

void NewtonMeshMatcherWrite::setup()
{
  const CFuint size = PhysicalModelStack::getActive()->getDim() - 1;
  _mappedCoord.resize(size);
  _mappedCoordMin.resize(size);
  _bestCoord.resize(size);
  _residual.resize(size);
  _otherResidual.resize(size);
  _jacobian.resize(size);
  // Reference Values for the Numerical Jacobian
  _refValues.resize(size);
  _refValues = 1.;

  _faceSelection.resize(_nbSelectedFaces);
  _tempCoord.resize(PhysicalModelStack::getActive()->getDim());
  _centroid.resize(PhysicalModelStack::getActive()->getDim());
}

//////////////////////////////////////////////////////////////////////////////

void NewtonMeshMatcherWrite::executeWrite(const CFuint iProc)
{

  CFAUTOTRACE;

  if(_shiftStates) shiftStatesCoord();


  // Names of the interfaces, trs
  const std::string interfaceName = getCommandGroupName();
  vector<std::string> otherTrsNames = getMethodData().getCoupledSubSystemsTRSNames(interfaceName);

  // If the geometry is non matching, then get the threshold for acceptance
  // else set the threshold to a very high value
  const bool isNonMatchingGeometry = getMethodData().isInterfaceGeometryNonMatching(interfaceName);
  CFreal nonMatchingGeometryThreshold = MathTools::MathConsts::CFrealMax();
  if(isNonMatchingGeometry) nonMatchingGeometryThreshold = getMethodData().getNonMatchingGeometryThreshold(interfaceName);

  // Get the geometric entity builder
  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  for (CFuint iTRS=0; iTRS < otherTrsNames.size(); iTRS++)
  {
    CFLogInfo("Matching for TRS [" << otherTrsNames[iTRS] << "]\n");
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

      cf_assert(interfaceCoords.size() == isAccepted.size());

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
        CFuint idx = iState - rejectedStates;
        bool accepted = false;
        _minimumDistanceOnFace = MathTools::MathConsts::CFrealMax();

        _coord.resize(PhysicalModelStack::getActive()->getDim());
        SubSysCouplerData::GeoEntityIdx matchingFace;

        const CFuint otherDim = interfaceCoords[iState].size();
        _coord = 0.;
        for(CFuint iDim = 0; iDim < otherDim; iDim++)
        {
          if(iDim < PhysicalModelStack::getActive()->getDim()) _coord[iDim] = interfaceCoords[iState][iDim];
        }

        // PreSelect a given number of faces whose centroid is close to the point to project
        facePreSelection();

        // Pair each fluid point with the closest wet structural element
        nodeToElementPairing(matchingFace);

        // Acceptation of the projection depending on
        // the distance between point and its projection
        if(nonMatchingGeometryThreshold > _minimumDistanceOnFace) accepted = true;

        if(accepted)
        {
          cf_assert(matchingFace.first.isNotNull());

          CFLogDebugMin("Accepted: " <<_coord << "\n");
          geoData.trs = matchingFace.first;
          geoData.idx = matchingFace.second;
          GeometricEntity& matchingFaceGeo = *geoBuilder->buildGE();

          (*coupledGeoEntities)[idx].first = matchingFace.first;
          (*coupledGeoEntities)[idx].second = matchingFace.second;

          // Compute the shape function at the projected point
          CFuint nbStatesInFace = matchingFaceGeo.getStates()->size();
          (*coupledGeoEntities)[idx].third.resize(nbStatesInFace);
          (*coupledGeoEntities)[idx].fourth.resize(PhysicalModelStack::getActive()->getDim());

          (*coupledGeoEntities)[idx].third = matchingFaceGeo.computeShapeFunctionAtMappedCoord(_mappedCoordMin);
          _tempCoord = matchingFaceGeo.computeCoordFromMappedCoord(_mappedCoordMin);
          (*coupledGeoEntities)[idx].fourth = _tempCoord;

          CFLogDebugMin("...and associated to: " << _tempCoord << "\n");
          CFLogDebugMin("...idx: " << idx << "\n");
          CFLogDebugMin("...Distance: " << _minimumDistanceOnFace << "\n");
          CFLogDebugMin("...ShapeFunction: " << (*coupledGeoEntities)[idx].third << "\n");
          std::vector<State*>* states = matchingFaceGeo.getStates();
          CFLogDebugMin("...nbStates: " << states->size() << "\n");
          CFLogDebugMin("...State[0]: " << (*states)[0]->getCoordinates() << "\n");
          CFLogDebugMin("...State[1]: " << (*states)[1]->getCoordinates() << "\n");
          if(states->size() > 2) { CFLogDebugMin("...State[2]: " << (*states)[2]->getCoordinates() << "\n"); }
          if(states->size() > 3) { CFLogDebugMin("...State[3]: " << (*states)[3]->getCoordinates() << "\n"); }

          //release the geometric entity
          geoBuilder->releaseGE();

          isAccepted[iState] = _minimumDistanceOnFace;

        }
        else
        {
          CFout << " =================== Rejected: " <<_coord << "\n";
          isAccepted[iState] = -1.;
          rejectedStates++;
        }
      } //end loop over otherStates

      //We know the number of accepted states -> can now resize the datahandle
      interfaceData.resize(otherNbStates - rejectedStates);

      //Write the file with the info isAccepted
      writeIsAcceptedFile(socketAcceptNames[iType]);
    } //end of loop over data transfer coord type
  } // end loop over the OtherTRS
  if(_shiftStates) unshiftStatesCoord();
}

//////////////////////////////////////////////////////////////////////////////

void NewtonMeshMatcherWrite::shiftStatesCoord()
{
  // Loop over all the nodes of all the TRS's of this command
  vector< SafePtr<TopologicalRegionSet> > trs = getTrsList();
  vector< Common::SafePtr<TopologicalRegionSet> >::iterator iTRS;

  /// @todo remove this this is very dirty!!!!!!!
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

void NewtonMeshMatcherWrite::unshiftStatesCoord()
{
  // Loop over all the nodes of all the TRS's of this command
  vector< SafePtr<TopologicalRegionSet> > trs = getTrsList();
  vector< Common::SafePtr<TopologicalRegionSet> >::iterator iTRS;

  /// @todo remove this this is very dirty!!!!!!!
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

void NewtonMeshMatcherWrite::facePreSelection()
{

  CFreal maxDistanceInSelection = MathTools::MathConsts::CFrealMax();
  CFreal distance;
  CFuint iMax = 0;

  //Initialize _faceSelection
  for(CFuint iFace=0;iFace < _nbSelectedFaces;iFace++)
  {
    _faceSelection[iFace].first = CFNULL;
  }

  // Loop over all the faces of all the TRS's of this command
  vector< SafePtr<TopologicalRegionSet> > trs = getTrsList();
  vector< Common::SafePtr<TopologicalRegionSet> >::iterator iTRS;

  //if we match against InnerCells, then resize all the temp Vectors
  for (iTRS = trs.begin(); iTRS != trs.end(); ++iTRS)
  {
    const std::string currentTRSName = (*iTRS)->getName();
    if(currentTRSName == "InnerCells"){
      cf_assert(trs.size() == 1);
      const CFuint size = PhysicalModelStack::getActive()->getDim();
      _mappedCoord.resize(size);
      _mappedCoordMin.resize(size);
      _bestCoord.resize(size);
      _residual.resize(size);
      _otherResidual.resize(size);
      _jacobian.resize(size);
      // Reference Values for the Numerical Jacobian
      _refValues.resize(size);
      _refValues = 1.;
    }
  }

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  for (iTRS = trs.begin(); iTRS != trs.end(); ++iTRS)
  {
    const CFuint nbGeos = (*iTRS)->getLocalNbGeoEnts();
    const std::string currentTRSName = (*iTRS)->getName();

    for(CFuint iGeoEnt = 0; iGeoEnt < nbGeos; ++iGeoEnt)
    {
      // build the GeometricEntity
      geoData.trs = (*iTRS);
      geoData.idx = iGeoEnt;

      _currentFace = &(*geoBuilder->buildGE());

      // Compute distance from point to centroid of the face
      _centroid = _currentFace->computeCentroid();

      distance = 0.;
      for(CFuint iDim = 0;iDim < _coord.size();iDim++)
      {
        distance += (_centroid[iDim]-_coord[iDim])*(_centroid[iDim]-_coord[iDim]);
      }

      //release the GeometricEntity
      geoBuilder->releaseGE();

      // If distance is smaller than the maximum of the preselected distance
      // then unselect the previous maximum and select this face
      // Then, recompute wich face is the further away...
      if (distance < maxDistanceInSelection)
      {
        //Set the new value at the position of the previous face
        _faceSelection[iMax].first = (*iTRS);
        _faceSelection[iMax].second = iGeoEnt;

        //Recompute the furthest face from the point
        maxDistanceInSelection = 0.;
        bool hasNull = false;
        for(CFuint iFace=0;iFace < _nbSelectedFaces;iFace++)
        {
          if(_faceSelection[iFace].first.isNotNull())
          {
            geoData.trs = _faceSelection[iFace].first;
            geoData.idx = _faceSelection[iFace].second;

            GeometricEntity& iFaceGeo = *geoBuilder->buildGE();

            _centroid = iFaceGeo.computeCentroid();
            distance = 0.;
            for(CFuint iDim = 0;iDim < _coord.size();iDim++)
            {
              distance += (_centroid[iDim]-_coord[iDim])*(_centroid[iDim]-_coord[iDim]);
            }

            if((distance>maxDistanceInSelection) && (!hasNull))
            {
              iMax = iFace;
              maxDistanceInSelection = distance;
            }

            // reset the trs to the current trs
            geoData.trs = (*iTRS);
            // release the geometric entity
            geoBuilder->releaseGE();
          }
          else
          {
            hasNull = true;
            iMax = iFace;
            maxDistanceInSelection = MathTools::MathConsts::CFrealMax();
          }
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void NewtonMeshMatcherWrite::nodeToElementPairing(SubSysCouplerData::GeoEntityIdx& matchingFace)
{

  Common::SafePtr<GeometricEntityPool<StdTrsGeoBuilder> >
  geoBuilder = getMethodData().getStdTrsGeoBuilder();

  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  // Loop over all the preSelected faces
  for (CFuint iFace = 0; iFace< _faceSelection.size(); ++iFace)
  {
//  CFout << "Face: " << iFace << "\n";
    if(_faceSelection[iFace].first.isNotNull())
    {
      geoData.trs = _faceSelection[iFace].first;
      geoData.idx = _faceSelection[iFace].second;

      _currentFace = &(*geoBuilder->buildGE());

//   CFout << "=====================================" << "\n";
//   CFout << "Node 0: " << *(_currentFace->getNode(0)) << "\n";
//   CFout << "Node 1: " << *(_currentFace->getNode(1)) << "\n";
//   CFout << "Node 2: " << *(_currentFace->getNode(2)) << "\n";
//   CFout << "=====================================" << "\n";
      // Check if node defined by coord has his projection in face...
      computeProjection();

      if (_bestDistance < _minimumDistanceOnFace)
      {
        matchingFace.first = geoData.trs;
        matchingFace.second = geoData.idx;

        _mappedCoordMin = _bestCoord;
        _minimumDistanceOnFace = _bestDistance;
// CFout << "Minimum Distance: " << _bestDistance << "\n";
      }
      // release the geometric entity
      geoBuilder->releaseGE();
    }
  } // end of loop over faces
}

///Old routine, without face preselection...
// void NewtonMeshMatcherWrite::nodeToElementPairing(GeometricEntity** matchingFace)
// {
//
// //  bool isOnCurrentFace(false);
// //  bool isOnAFace(false);
//
//   /// Loop over all the faces of all the TRS's of this command
//   vector<TopologicalRegionSet*>* trs = getTrsList();
//   vector< Common::SafePtr<TopologicalRegionSet> >::iterator iTRS;
//   for (iTRS = trs.begin(); iTRS != trs.end(); ++iTRS) {
//     const CFuint nbTRs = (*iTRS)->getNbTRs();
//     // loop over all the topological regions
//     for (CFuint iTR = 0; iTR < nbTRs; ++iTR) {
//       GeomEntList *const faces = (*(*iTRS))[iTR]->getGeomEntList();
//       GeomEntList::iterator itg;
//       for (itg = faces->begin(); itg != faces->end(); ++itg) {
//         _currentFace = *itg;
//
//         /// Check if node defined by coord has his projection in face...
//         computeProjection();
//
//         if (_residual[0] < _minimumDistanceOnFace)
//         {
//           *matchingFace = _currentFace;
//           _mappedCoordMin.resize(_mappedCoord.size());
//           _mappedCoordMin = _mappedCoord;
//           _minimumDistanceOnFace = _residual[0];
//         }
//       } // end of loop over faces
//     } // end of loop over TR's
//   } // end of loop over TRS's
// }


//////////////////////////////////////////////////////////////////////////////

void NewtonMeshMatcherWrite::computeProjection()
{
  setupNewton();

  bool isAchieved = false;
  for(CFuint k = 0; !isAchieved; ++k)
  {
    //Compute jacobian
    takeStep();

    //Update the solution
    newUpdateSolution();

    //compute the new distance
    const CFreal newDistance = computeResidual();

    if(newDistance < _bestDistance)
    {
      _bestCoord = _mappedCoord;
      _bestDistance = newDistance;
// CFout << "***********************************\n";
// CFout << "_bestCoord: " << _bestCoord <<"\n";
// CFout << "_bestDistance: " << _bestDistance <<"\n";
    }

    //stop condition: difference between two iterations and absolute number of steps
    if((k >= _maxNewtonSteps) || (_bestDistance< _maxAccuracy))
    {
      isAchieved = true;
// CFout << "*********FINISHED**************\n";
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void NewtonMeshMatcherWrite::setupNewton()
{
  // Setup
  // We have to minimize the distance from the point to the face

  _numericalJacob->setRefValues(_refValues);

  // Initialization
  _mappedCoord = 0.25;

  _currentIVar = 0;

  _bestCoord = 0.;
  _bestDistance = MathTools::MathConsts::CFrealMax();
}

//////////////////////////////////////////////////////////////////////////////

CFreal NewtonMeshMatcherWrite::computeResidual()
{
  _tempCoord = _currentFace->computeCoordFromMappedCoord(_mappedCoord);

  CFreal d2 = 0;
  for(CFuint i=0;i<_tempCoord.size();i++)
  {
    d2 += (_coord[i]-_tempCoord[i])*(_coord[i]-_tempCoord[i]);
  }

  return sqrt(d2);
}

//////////////////////////////////////////////////////////////////////////////

void NewtonMeshMatcherWrite::takeStep()
{

  for (CFuint iVar = 0; iVar < _mappedCoord.size(); ++iVar) {
    _residual[iVar] = computeResidual();
    _numericalJacob->perturb(iVar, _mappedCoord[iVar]);
    _otherResidual[iVar] = computeResidual();
    _numericalJacob->restore(_mappedCoord[iVar]);
  }

  // jacobian contribution (dR/dU)_k
  // compute (R[U + dU] - R[U])/eps
  _numericalJacob->computeDerivative(_residual,
                                     _otherResidual,
                                     _jacobian);
}

//////////////////////////////////////////////////////////////////////////////

void NewtonMeshMatcherWrite::newUpdateSolution()
{

//  std::cout << "Jacob: " << _jacobian<< std::endl;
  CFreal maxGradient = -MathTools::MathConsts::CFrealMax();
  for (CFuint iVar = 0; iVar < _jacobian.size(); ++iVar)
  {
    if(fabs(_jacobian[iVar]) > maxGradient){
      _currentIVar = iVar;
      maxGradient = fabs(_jacobian[iVar]);
    }
  }

  //At this step, we update _currentIVar component
  const CFuint nbNodes = _currentFace->getNbNodesGeometryShapeFunction();

  // Define the constraints
  CFreal maxi;
  CFreal mini;

  if(nbNodes==3) //triangle
  {
    cf_assert(_mappedCoord.size() == 2);
    mini = 0.;
    if (_currentIVar == 0) maxi = 1.-_mappedCoord[1];
    else maxi = 1.-_mappedCoord[0];
  }
  else //line or quad
  {
    maxi = 1.;
    mini = -1.;
  }

  CFreal tolerance = 1.e-10;
  if(((_mappedCoord[_currentIVar] - mini) < tolerance) && (_jacobian[_currentIVar] > 0.)){
    const CFuint nbVars = _mappedCoord.size();
    if(_currentIVar == nbVars-1) _currentIVar = 0;
    else _currentIVar++;
  }
  else{
    if(((maxi - _mappedCoord[_currentIVar]) < tolerance) && (_jacobian[_currentIVar] < 0.)){
      const CFuint nbVars = _mappedCoord.size();
      if(_currentIVar == nbVars-1) _currentIVar = 0;
      else _currentIVar++;
    }
  }

  //while updated value is not feasible do
  //  if updated value is within the feasible region then update
  //  else then modify the step by half
  bool isFeasible(false);
  CFreal dU = -_residual[_currentIVar] / _jacobian[_currentIVar];
  CFreal originalDistance = computeResidual();
  CFreal relax = 0.5;

  while(!isFeasible)
  {
    CFreal last = _mappedCoord[_currentIVar];
    _mappedCoord[_currentIVar] += relax*dU;
    if ((_mappedCoord[_currentIVar] > (maxi+MathTools::MathConsts::CFrealEps())) ||
        (_mappedCoord[_currentIVar] < (mini-MathTools::MathConsts::CFrealEps())) ||
        (computeResidual() > originalDistance))
    {
      relax *= 0.5;
      if(std::fabs(_mappedCoord[_currentIVar] - last) < _maxAccuracy) isFeasible = true;
      _mappedCoord[_currentIVar] = last;
    }
    else
    {
      isFeasible = true;
    }
  }

// CFout <<"_mappedCoord: " << _mappedCoord <<"\n";
  //Change the update variable to iVar+1
  // if (iVar+1) > nbVars then set iVar = 0
//   const CFuint nbVars = _mappedCoord.size();
//   if(_currentIVar == nbVars-1) _currentIVar = 0;
//   else _currentIVar++;

//std::cout << "Distance: " << computeResidual() << " - " << _mappedCoord[0] << " - " << _mappedCoord[1] << std::endl;

}

//////////////////////////////////////////////////////////////////////////////

///old way of updating
void NewtonMeshMatcherWrite::updateSolution()
{
  CFreal dU;
  CFreal maxi;
  CFuint nbNodes = _currentFace->getNbNodesSolutionShapeFunction();

  for(CFuint iVar=0;iVar<_mappedCoord.size();++iVar)
  {
    dU = -_residual[iVar] / _jacobian[iVar];
    _mappedCoord[iVar] += 0.1*dU;

    //Add constraints
    if(nbNodes==3) //triangle
    {
      cf_assert(_mappedCoord.size() == 2);
      if (iVar == 0) maxi = 1.-min(1.,max(0.,_mappedCoord[0]));
      else maxi = 1.-min(1.,max(0.,_mappedCoord[1]));

      _mappedCoord[iVar] = min(maxi,max(0.,_mappedCoord[iVar]));
    }
    else //line or quad
    {
      _mappedCoord[iVar] = min(1.,max(-1.,_mappedCoord[iVar]));
    }

  }

  for(CFuint iVar=0;iVar<_mappedCoord.size();++iVar)
  {
    _residual[iVar] = computeResidual();
  }

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshMatcher

  } // namespace Numerics

} // namespace COOLFluiD
