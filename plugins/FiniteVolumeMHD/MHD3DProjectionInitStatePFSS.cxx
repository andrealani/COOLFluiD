#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"

#include "FiniteVolumeMHD/MHD3DProjectionInitStatePFSS.hh"
#include "Common/CFLog.hh"
#include "Common/PEFunctions.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/DataHandle.hh"
#include "Common/BadValueException.hh"
#include "Common/FilesystemException.hh"
#include "MHD/MHD3DProjectionVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<MHD3DProjectionInitStatePFSS, CellCenterFVMData, FiniteVolumeMHDModule>
mHD3DProjectionInitStatePFSSProvider("MHD3DProjectionInitStatePFSS");

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionInitStatePFSS::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("epsilon","Desired accuracy for Parker solution.");
   options.addConfigOption< CFreal >("rMin","Radius of the inner boundary.");
   options.addConfigOption< CFreal >("rMax","Radius of the outer boundary.");
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionInitStatePFSS::MHD3DProjectionInitStatePFSS(const std::string& name) :
  InitState(name),
  _varSet(CFNULL),
  socket_AlmrealBegin("AlmrealBegin"),
  socket_AlmimgBegin("AlmimgBegin"),
  socket_BlmrealBegin("BlmrealBegin"),
  socket_BlmimgBegin("BlmimgBegin"),
  socket_AlmrealEnd("AlmrealEnd"),
  socket_AlmimgEnd("AlmimgEnd"),
  socket_BlmrealEnd("BlmrealEnd"),
  socket_BlmimgEnd("BlmimgEnd"),
  socket_BPFSS("BPFSS"),
  _nameBeginPFSSDataFile(),
  _nameEndPFSSDataFile()
{
  addConfigOptionsTo(this);

  _epsilon = 1.0e-6;
  setParameter("epsilon",&_epsilon);

  _rMin = 1.03;
  setParameter("rMin",&_rMin);

  _rMax = 20.0;
  setParameter("rMax",&_rMax);

}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionInitStatePFSS::~MHD3DProjectionInitStatePFSS()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionInitStatePFSS::readInputFile()
{
  DataHandle<std::vector<CFreal> > AlmrealBegin  = socket_AlmrealBegin.getDataHandle();
  DataHandle<std::vector<CFreal> > AlmimgBegin  = socket_AlmimgBegin.getDataHandle();
  DataHandle<std::vector<CFreal> > BlmrealBegin  = socket_BlmrealBegin.getDataHandle();
  DataHandle<std::vector<CFreal> > BlmimgBegin  = socket_BlmimgBegin.getDataHandle();
  DataHandle<std::vector<CFreal> > AlmrealEnd  = socket_AlmrealEnd.getDataHandle();
  DataHandle<std::vector<CFreal> > AlmimgEnd  = socket_AlmimgEnd.getDataHandle();
  DataHandle<std::vector<CFreal> > BlmrealEnd  = socket_BlmrealEnd.getDataHandle();
  DataHandle<std::vector<CFreal> > BlmimgEnd  = socket_BlmimgEnd.getDataHandle();

  // The order of coefficients should be Almreal, Almimg, Blmreal and Blmimg each separated with a blank line in the input file 

  // First file

  SelfRegistPtr<Environment::FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  ifstream& inputFile = fhandle->open(constructFilename(_nameBeginPFSSDataFile));

  std::string line, tmpstr;
  vector<CFreal> tmprow;

  const CFuint nbLModes = _varSet->getNbLModes();

  for (CFuint iLine = 0; iLine <= nbLModes; ++iLine) {
      getline(inputFile,line);
      istringstream liness( line );
  
      for (CFuint jLine = 0; jLine <= nbLModes; ++jLine) {
       	getline(liness, tmpstr, ' ');
       	tmprow.push_back(atof(tmpstr.c_str())); 
      }
      AlmrealBegin[iLine] = tmprow;
      tmprow.clear();
  }

  // skip the undesired modes and the blank line separating the coefficients
  while (getline(inputFile,line)) {
        if (line == "") {
           break;
        }
  }

  for (CFuint iLine = 0; iLine <= nbLModes; ++iLine) {
      getline(inputFile,line);
      istringstream liness( line );

      for (CFuint jLine = 0; jLine <= nbLModes; ++jLine) {
        getline(liness, tmpstr, ' ');
        tmprow.push_back(atof(tmpstr.c_str()));
      }
      AlmimgBegin[iLine] = tmprow;
      tmprow.clear();
  }

  // skip the undesired modes and the blank line separating the coefficients
  while (getline(inputFile,line)) {
        if (line == "") {
           break;
        }
  }

  for (CFuint iLine = 0; iLine <= nbLModes; ++iLine) {
      getline(inputFile,line);
      istringstream liness( line );

      for (CFuint jLine = 0; jLine <= nbLModes; ++jLine) {
        getline(liness, tmpstr, ' ');
        tmprow.push_back(atof(tmpstr.c_str()));
      }
      BlmrealBegin[iLine] = tmprow;
      tmprow.clear();
  }

  // skip the undesired modes and the blank line separating the coefficients
  while (getline(inputFile,line)) {
        if (line == "") {
           break;
        }
  }

  for (CFuint iLine = 0; iLine <= nbLModes; ++iLine) {
      getline(inputFile,line);
      istringstream liness( line );

      for (CFuint jLine = 0; jLine <= nbLModes; ++jLine) {
        getline(liness, tmpstr, ' ');
        tmprow.push_back(atof(tmpstr.c_str()));
      }
      BlmimgBegin[iLine] = tmprow;
      tmprow.clear();
  }

  fhandle->close();

  // Second file

  ifstream& inpFile = fhandle->open(constructFilename(_nameEndPFSSDataFile));

  for (CFuint iLine = 0; iLine <= nbLModes; ++iLine) {
      getline(inpFile,line);
      istringstream liness( line );

      for (CFuint jLine = 0; jLine <= nbLModes; ++jLine) {
        getline(liness, tmpstr, ' ');
        tmprow.push_back(atof(tmpstr.c_str()));
      }
      AlmrealEnd[iLine] = tmprow;
      tmprow.clear();
  }

  // skip the undesired modes and the blank line separating the coefficients
  while (getline(inpFile,line)) {
        if (line == "") {
           break;
        }
  }

  for (CFuint iLine = 0; iLine <= nbLModes; ++iLine) {
      getline(inpFile,line);
      istringstream liness( line );

      for (CFuint jLine = 0; jLine <= nbLModes; ++jLine) {
        getline(liness, tmpstr, ' ');
        tmprow.push_back(atof(tmpstr.c_str()));
      }
      AlmimgEnd[iLine] = tmprow;
      tmprow.clear();
  }

  // skip the undesired modes and the blank line separating the coefficients
  while (getline(inpFile,line)) {
        if (line == "") {
           break;
        }
  }

  for (CFuint iLine = 0; iLine <= nbLModes; ++iLine) {
      getline(inpFile,line);
      istringstream liness( line );

      for (CFuint jLine = 0; jLine <= nbLModes; ++jLine) {
        getline(liness, tmpstr, ' ');
        tmprow.push_back(atof(tmpstr.c_str()));
      }
      BlmrealEnd[iLine] = tmprow;
      tmprow.clear();
  }

  // skip the undesired modes and the blank line separating the coefficients
  while (getline(inpFile,line)) {
        if (line == "") {
           break;
        }
  }

  for (CFuint iLine = 0; iLine <= nbLModes; ++iLine) {
      getline(inpFile,line);
      istringstream liness( line );

      for (CFuint jLine = 0; jLine <= nbLModes; ++jLine) {
        getline(liness, tmpstr, ' ');
        tmprow.push_back(atof(tmpstr.c_str()));
      }
      BlmimgEnd[iLine] = tmprow;
      tmprow.clear();
  }

  fhandle->close();

}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionInitStatePFSS::computeParkerSolution(const CFreal r, 
					  RealVector& velParkerSpherical,
					  CFreal& rhoParker,
					  CFreal& pParker)
{
  // gravitational constant
  const CFreal G = 6.67384e-11;

  // Boltzmann constant
  const CFreal k = 1.3806503e-23;

  // mass of proton
  const CFreal mp = 1.67262158e-27;

  const CFreal TRef = _varSet->getTRef();
  const CFreal lRef = _varSet->getLRef();

  const CFreal mass = _varSet->getMass();

  const CFreal vRef = sqrt(2.0*k*TRef/mp);
 
  // calculate the escape velocity: sqrt(2GM/R)
  const CFreal vEsc = sqrt(2.0*G*mass/lRef)/vRef;

  // calculate the isothermal speed of sound
  const CFreal a = sqrt(2.0*k*TRef/mp)/vRef;

  // calculate the expected transsonic point
  const CFreal rTrans = (vEsc*vEsc)/(4.0*a*a);

  // check whether the transsonic point is in the computational domain
  cf_assert(rTrans > _rMin);
  cf_assert(rTrans < _rMax);
  cf_assert(r != rTrans);

  // Parker's solar wind is radial
  velParkerSpherical[1] = 0.0;
  velParkerSpherical[2] = 0.0; 

  CFreal ur0 = 1.0, ur1, urh, diff;
  bool hasSolution = false;

  // calculate the Parker's solar wind speed that goes from subsonic to supersonic regime
  if (r < rTrans) {
	for (CFuint iter = 1; iter < 1000; ++iter) {	
	    //ur1 = ((vEsc*vEsc*vEsc*vEsc)/(16.0*r*r))*exp(0.5*((ur0*ur0)+3.0-((vEsc*vEsc)/r)));
	    ur1 = ((rTrans*rTrans)/(r*r))*exp(0.5*((ur0*ur0)+3.0-(4.0*rTrans/r)));
            diff = fabs(ur1-ur0);
	    if (diff < _epsilon) {
		velParkerSpherical[0] = ur1*a;
                hasSolution = true;
                break;
            }
	    else {
	        ur0 = ur1;
	    }
     	}
        cf_assert(hasSolution);
  }
  if (r > rTrans) {
        for (CFuint iter = 1; iter < 1000; ++iter) {
            //urh = ((vEsc*vEsc)/r)-3.0+2.0*log(16.0*ur0*r*r/(vEsc*vEsc*vEsc*vEsc));
	    urh = (4.0*rTrans/r)-3.0+2.0*log(ur0*r*r/(rTrans*rTrans));
	    ur1 = sqrt(urh);
            diff = fabs(ur1-ur0);
            if (diff < _epsilon) {
                velParkerSpherical[0] = ur1*a;
                hasSolution = true;
                break;
            }
            else {
                ur0 = ur1;
            }
        }
        cf_assert(hasSolution);
  }
  
  // calculate the approximate base solar wind speed in the vicinity of the photosphere (i.e. at _rMin)
  const CFreal vBase = a*((rTrans*rTrans)/(_rMin*_rMin))*exp(1.5-(2.0*rTrans/_rMin));
  // density obeys rho*v_r*r^2 = const.
  // assuming rhoBase = 1
  rhoParker = (vBase*_rMin*_rMin)/(velParkerSpherical[0]*r*r);

  // relation between density and pressure is p=a*a*rho
  pParker = a*a*rhoParker;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionInitStatePFSS::executeOnTrs()
{
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "InitState::executeOnTrs() called for TRS: "
  << trs->getName() << "\n");

  if (trs->getName() != "InnerFaces") {
    throw BadValueException (FromHere(),"InitState not applied to InnerFaces!!!");
  }

  // this cannot be used for FV boundary faces because
  // ghost state and inner state could have the same local ID !!!
  SafePtr<vector<CFuint> > trsStates = trs->getStatesInTrs();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  std::vector<CFuint>::iterator itd;
  CFreal rhoParker,pParker;
  RealVector stateCoordsCartesian(PhysicalModelStack::getActive()->getDim()),
	stateCoordsSpherical(PhysicalModelStack::getActive()->getDim()), 
	BPFSSCartesian(PhysicalModelStack::getActive()->getDim()), 
	velParkerSpherical(PhysicalModelStack::getActive()->getDim()), 
	velParkerCartesian(PhysicalModelStack::getActive()->getDim());
  RealMatrix carSphTransMat(PhysicalModelStack::getActive()->getDim(),PhysicalModelStack::getActive()->getDim()), 
	sphCarTransMat(PhysicalModelStack::getActive()->getDim(),PhysicalModelStack::getActive()->getDim());
  if(_inputAdimensionalValues)
  {
    for (itd = trsStates->begin(); itd != trsStates->end(); ++itd) {
      State* const currState = states[(*itd)];
      stateCoordsCartesian = currState->getCoordinates();
      _varSet->setTransformationMatrices(stateCoordsCartesian,stateCoordsSpherical,carSphTransMat,sphCarTransMat); 
      computeParkerSolution(stateCoordsSpherical[0],velParkerSpherical,rhoParker,pParker);
      velParkerCartesian = sphCarTransMat*velParkerSpherical;

      const CFreal gammaMinus1 = _varSet->getModel()->getGamma() - 1.0;

      const CFreal rhoE1Parker = (pParker/gammaMinus1) + 0.5*rhoParker*(velParkerCartesian[0]*velParkerCartesian[0]+
									velParkerCartesian[1]*velParkerCartesian[1]+
									velParkerCartesian[2]*velParkerCartesian[2]);

      // input variables are conservative
      (*_input)[0] = rhoParker;
      (*_input)[1] = rhoParker*velParkerCartesian[0];
      (*_input)[2] = rhoParker*velParkerCartesian[1];
      (*_input)[3] = rhoParker*velParkerCartesian[2];
      //!!! VARIABLE MAGNETIC FIELD B1 IS ASSIGNED TO BE ZERO SINCE THE INITIAL TOTAL MAGNETIC FIELD ENTIRELY CONSISTS OF THE POTENTIAL FIELD B0
      (*_input)[4] = 0.0;
      (*_input)[5] = 0.0;
      (*_input)[6] = 0.0;
      (*_input)[7] = rhoE1Parker;
      (*_input)[8] = 0.0;
      *currState = *_inputToUpdateVar->transform(_input);
    }
  }
  else
  {
    State dimState;
    for (itd = trsStates->begin(); itd != trsStates->end(); ++itd) {
      State* const currState = states[(*itd)];
      stateCoordsCartesian = currState->getCoordinates();
      _varSet->setTransformationMatrices(stateCoordsCartesian,stateCoordsSpherical,carSphTransMat,sphCarTransMat);
      computeParkerSolution(stateCoordsSpherical[0],velParkerSpherical,rhoParker,pParker);
      velParkerCartesian = sphCarTransMat*velParkerSpherical;

      const CFreal gammaMinus1 = _varSet->getModel()->getGamma() - 1.0;
      
      const CFreal rhoE1Parker = (pParker/gammaMinus1) + 0.5*rhoParker*(velParkerCartesian[0]*velParkerCartesian[0]+
                                                                        velParkerCartesian[1]*velParkerCartesian[1]+
                                                                        velParkerCartesian[2]*velParkerCartesian[2]);

      // input variables are conservative
      (*_input)[0] = rhoParker;
      (*_input)[1] = rhoParker*velParkerCartesian[0];
      (*_input)[2] = rhoParker*velParkerCartesian[1];
      (*_input)[3] = rhoParker*velParkerCartesian[2];
      //!!! VARIABLE MAGNETIC FIELD B1 IS ASSIGNED TO BE ZERO SINCE THE INITIAL TOTAL MAGNETIC FIELD ENTIRELY CONSISTS OF THE POTENTIAL FIELD B0
      (*_input)[4] = 0.0;
      (*_input)[5] = 0.0;
      (*_input)[6] = 0.0;
      (*_input)[7] = rhoE1Parker;
      (*_input)[8] = 0.0;
      dimState = *_inputToUpdateVar->transform(_input);
      _varSet->setAdimensionalValues(dimState, *currState);
    }
  }
}


//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionInitStatePFSS::setup()
{

  InitState::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DProjectionVarSet>();
  cf_assert(_varSet.isNotNull());

  const std::string potentialBType = _varSet->getPotentialBType();

  if (potentialBType == "PFSS") {

	DataHandle<std::vector<CFreal> > AlmrealBegin  = socket_AlmrealBegin.getDataHandle();
        DataHandle<std::vector<CFreal> > AlmimgBegin  = socket_AlmimgBegin.getDataHandle();
        AlmrealBegin.resize(_varSet->getNbLModes()+1);
        AlmimgBegin.resize(_varSet->getNbLModes()+1);

        DataHandle<std::vector<CFreal> > BlmrealBegin  = socket_BlmrealBegin.getDataHandle();
        DataHandle<std::vector<CFreal> > BlmimgBegin  = socket_BlmimgBegin.getDataHandle();
        BlmrealBegin.resize(_varSet->getNbLModes()+1);
        BlmimgBegin.resize(_varSet->getNbLModes()+1);

        DataHandle<std::vector<CFreal> > AlmrealEnd  = socket_AlmrealEnd.getDataHandle();
        DataHandle<std::vector<CFreal> > AlmimgEnd  = socket_AlmimgEnd.getDataHandle();
        AlmrealEnd.resize(_varSet->getNbLModes()+1);
        AlmimgEnd.resize(_varSet->getNbLModes()+1);

        DataHandle<std::vector<CFreal> > BlmrealEnd  = socket_BlmrealEnd.getDataHandle();
        DataHandle<std::vector<CFreal> > BlmimgEnd  = socket_BlmimgEnd.getDataHandle();
        BlmrealEnd.resize(_varSet->getNbLModes()+1);
        BlmimgEnd.resize(_varSet->getNbLModes()+1);

	bool interpolationFlag = false;
        _varSet->getModel()->setInterpolationFlag(interpolationFlag);

        _nameBeginPFSSDataFile = _varSet->getNameBeginPFSSCoeffFile();
        _nameEndPFSSDataFile = _varSet->getNameEndPFSSCoeffFile();
   
	// if this is a parallel simulation, only ONE process at a time reads the file
  	const std::string nsp = this->getMethodData().getNamespace();
	runSerial<void, MHD3DProjectionInitStatePFSS, 
		  &MHD3DProjectionInitStatePFSS::readInputFile>(this, nsp);
	
	SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
        const CFuint nbTrsFaces = trs->getLocalNbGeoEnts();

        CFout << "PFSS magnetic field is computed for " << trs->getName() << " on " << nbTrsFaces << " faces\n";
	
        DataHandle<std::vector<CFreal> > BPFSS = socket_BPFSS.getDataHandle();
        BPFSS.resize(MeshDataStack::getActive()->Statistics().getNbFaces());

        Common::SafePtr<GeometricEntityPool<FaceTrsGeoBuilder> >
                geoBuilder = getMethodData().getFaceTrsGeoBuilder();

        SafePtr<FaceTrsGeoBuilder> geoBuilderPtr = geoBuilder->getGeoBuilder();
        geoBuilderPtr->setDataSockets(socket_states, socket_gstates,
                                socket_nodes);

        FaceTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();
        geoData.isBFace = false;
        geoData.trs = trs;

	CFreal quadPointXCoord, quadPointYCoord, quadPointZCoord;
        // DataHandle does not accept RealVector so temporary vectors are created
        RealVector quadPointCoordsCartesian(PhysicalModelStack::getActive()->getDim()),
                quadPointCoordsSpherical(PhysicalModelStack::getActive()->getDim()),
                BPFSSCartesian(PhysicalModelStack::getActive()->getDim());
        RealMatrix carSphTransMat(PhysicalModelStack::getActive()->getDim(),PhysicalModelStack::getActive()->getDim()),
                sphCarTransMat(PhysicalModelStack::getActive()->getDim(),PhysicalModelStack::getActive()->getDim());

        vector<CFreal> BPFSSCartesianCoords(PhysicalModelStack::getActive()->getDim());

	for (CFuint iFace = 0; iFace < nbTrsFaces; ++iFace) {
                CFLogDebugMed( "iFace = " << iFace << "\n");

                // build the GeometricEntity
                geoData.idx = iFace;

                GeometricEntity *const face = geoBuilder->buildGE();
                const vector<Node*>& nodes = *face->getNodes();

                const CFuint nbNodesInFace = nodes.size();

                quadPointXCoord = 0.0;
                quadPointYCoord = 0.0;
                quadPointZCoord = 0.0;

		for (CFuint iNode = 0; iNode < nbNodesInFace; ++iNode) {
                        quadPointXCoord += (*(nodes[iNode]))[XX];
                        quadPointYCoord += (*(nodes[iNode]))[YY];
                        quadPointZCoord += (*(nodes[iNode]))[ZZ];
                }

                quadPointXCoord /= nbNodesInFace;
                quadPointYCoord /= nbNodesInFace;
                quadPointZCoord /= nbNodesInFace;

                quadPointCoordsCartesian[0] = quadPointXCoord;
                quadPointCoordsCartesian[1] = quadPointYCoord;
                quadPointCoordsCartesian[2] = quadPointZCoord;

                const CFuint faceID = face->getID();

		_varSet->setTransformationMatrices(quadPointCoordsCartesian,quadPointCoordsSpherical,carSphTransMat,sphCarTransMat);
                _varSet->computePFSSMagneticField(quadPointCoordsSpherical,BPFSSCartesian,sphCarTransMat);

                BPFSSCartesianCoords[0] = BPFSSCartesian[0];
                BPFSSCartesianCoords[1] = BPFSSCartesian[1];
                BPFSSCartesianCoords[2] = BPFSSCartesian[2];

                BPFSS[faceID] = BPFSSCartesianCoords;

                // release the GeometricEntity
                geoBuilder->releaseGE();
        }

  }
}

//////////////////////////////////////////////////////////////////////////////

boost::filesystem::path MHD3DProjectionInitStatePFSS::constructFilename(std::string fileName)
{
  boost::filesystem::path fpath(fileName);

  CFout << "Reading PFSS spherical harmonics coefficients from: " << fpath.string() << "\n";

  return fpath;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
MHD3DProjectionInitStatePFSS::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result =
          InitState::providesSockets();
  result.push_back(&socket_AlmrealBegin);
  result.push_back(&socket_AlmimgBegin);
  result.push_back(&socket_BlmrealBegin);
  result.push_back(&socket_BlmimgBegin);
  result.push_back(&socket_AlmrealEnd);
  result.push_back(&socket_AlmimgEnd);
  result.push_back(&socket_BlmrealEnd);
  result.push_back(&socket_BlmimgEnd);
  result.push_back(&socket_BPFSS);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
