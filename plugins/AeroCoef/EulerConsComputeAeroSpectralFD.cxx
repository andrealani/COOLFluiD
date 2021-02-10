#include "AeroCoef/AeroCoefSpectralFD.hh"
#include "AeroCoef/EulerConsComputeAeroSpectralFD.hh"

#include "Common/PE.hh"

#include "Environment/DirPaths.hh"

#include "Common/BadValueException.hh"
#include "Framework/FaceToCellGEBuilder.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/PathAppender.hh"
#include "Common/ParserException.hh"
#include "Framework/SubSystemStatus.hh"

#include "NavierStokes/EulerVarSet.hh"

#include "SpectralFD/BaseBndFaceTermComputer.hh"
#include "SpectralFD/BCStateComputer.hh"
#include "SpectralFD/SpectralFDElementData.hh"
#include "SpectralFD/SpectralFDMethod.hh"
#include "SpectralFD/SpectralFDMethodData.hh"
#include "SpectralFD/TensorProductGaussIntegrator.hh"

#include "Common/SafePtr.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::SpectralFD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<EulerConsComputeAeroSpectralFD, DataProcessingData, AeroCoefSpectralFDModule>
    EulerConsComputeAeroSpectralFDProvider("EulerConsComputeAeroSpectralFD");

//////////////////////////////////////////////////////////////////////////////

void EulerConsComputeAeroSpectralFD::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("Alpha","Definition of the function defining the angle between the flow and the airfoil.");
   options.addConfigOption< CFuint >("SaveRateAero","Rate for saving the output file with aerodynamic coeficients.");
   options.addConfigOption< std::string >("OutputFileAero","Name of Output File to write the results.");
   options.addConfigOption< CFuint >("SaveRateWall","Save Output File containing the wall values every...iterations.");
   options.addConfigOption< std::string >("OutputFileWall","Name of Output File to write the wall values.");
   options.addConfigOption< CFreal >("uInf","Velocity at infinity.");
   options.addConfigOption< CFreal >("rhoInf","Density at infinity.");
   options.addConfigOption< CFreal >("pInf","Pressure at infinity.");
   options.addConfigOption< bool >("AppendTime","Append time to file name.");
   options.addConfigOption< bool >("AppendIter","Append Iteration# to file name.");

}

//////////////////////////////////////////////////////////////////////////////

EulerConsComputeAeroSpectralFD::EulerConsComputeAeroSpectralFD(const std::string& name) :
  DataProcessingCom(name),
  m_sockets(),
  m_faceBuilder(CFNULL),
  m_bndFaceTermComputer(CFNULL),
  m_bcStateComputer(CFNULL),
  m_updateVarSet(CFNULL),
  m_vars("t"),
  m_eval(0.0,1),
  m_lift(),
  m_drag(),
  m_alphadeg(0.),
  m_alpharad(0.)
{
  m_fileWall = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();

  addConfigOptionsTo(this);

  m_nameOutputFileWall = "Wall";
  setParameter("OutputFileWall",&m_nameOutputFileWall);

  m_saveRateWall = 1;
  setParameter("SaveRateWall",&m_saveRateWall);

  m_function = std::string("0.");
  setParameter("Alpha",&m_function);

  m_nameOutputFileAero = "AeroCoef.plt";
  setParameter("OutputFileAero",&m_nameOutputFileAero);

  m_saveRateAero = 1;
  setParameter("SaveRateAero",&m_saveRateAero);

  m_uInf = 0.;
  setParameter("uInf",&m_uInf);

  m_rhoInf = 0.;
  setParameter("rhoInf",&m_rhoInf);

  m_pInf = 0.;
  setParameter("pInf",&m_pInf);

  m_appendIter = false;
  setParameter("AppendIter",&m_appendIter);

  m_appendTime = false;
  setParameter("AppendTime",&m_appendTime);
}

//////////////////////////////////////////////////////////////////////////////

EulerConsComputeAeroSpectralFD::~EulerConsComputeAeroSpectralFD()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
EulerConsComputeAeroSpectralFD::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result = m_sockets.getAllSinkSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
EulerConsComputeAeroSpectralFD::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = m_sockets.getAllSourceSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void EulerConsComputeAeroSpectralFD::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);

  // Create the socket source for the boundaryNormals for each TRS to which this command applies
  // This is because we have different normals in FiniteVolume and FluctSplit / FiniteElement
//   const CFuint nbTrs = getTrsNames().size();
//   for (CFuint iTrs = 0; iTrs < nbTrs; ++iTrs)
//   {
//     std::string socketName = getTrsName(iTrs) + "-boundaryNormals";
//     m_sockets.createSocketSource<const CFreal*>(socketName);
//   }

  try {
    setFunction();
  }
  catch (Common::Exception& e) {
    CFout << e.what() << "\n";
    throw;
  }
}

//////////////////////////////////////////////////////////////////////////////

void EulerConsComputeAeroSpectralFD::setFunction()
{
  // some sanity checks
  std::vector<std::string> m_functionDef = Common::StringOps::getWords(m_function);
  if(m_functionDef.size() != 1)
   throw BadValueException (FromHere(),"EulerConsComputeAeroSpectralFD::setFuntion(): alpha function wrongly defined.");

  m_functionParser.Parse(m_function, m_vars);

  if (m_functionParser.ErrorMsg() != 0) {
    std::string msg("ParseError in CFL::setFuntion(): ");
    msg += std::string(m_functionParser.ErrorMsg());
    msg += " Function: " + m_function;
    msg += " Vars: "     + m_vars;
    throw Common::ParserException (FromHere(),msg);
  }
}

//////////////////////////////////////////////////////////////////////////////

void EulerConsComputeAeroSpectralFD::executeOnTrs()
{
  CFAUTOTRACE;

  CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  // compute the value of the angle
  m_eval[0] = SubSystemStatusStack::getActive()->getCurrentTime();
  m_alphadeg = m_functionParser.Eval(&m_eval[0]);

  // transform into Radiants
  m_alpharad = m_alphadeg*MathTools::MathConsts::CFrealPi()/180;

  // Execute and save file if needed...
  if(!(iter % m_saveRateWall))
  {
      computeWall();
  }

  if(!(iter % m_saveRateAero))
  {
    computeAero();
  }
}

//////////////////////////////////////////////////////////////////////////////

void EulerConsComputeAeroSpectralFD::computeFaceData()
{
  // set the face in the boundary face term computer
  m_bndFaceTermComputer->setCurrentFace(m_face);

  // compute the face data in the boundary face term computer
  m_bndFaceTermComputer->computeFacePntSetData();
}

//////////////////////////////////////////////////////////////////////////////

void EulerConsComputeAeroSpectralFD::computeWall()
{
  CFAUTOTRACE;

  prepareOutputFileWall(); // file handle is opened here

  // get some data needed further
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqsM1 = nbEqs - 1;
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFreal R = m_updateVarSet->getModel()->getR();
  const CFreal gamma = m_updateVarSet->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;

  cf_assert(m_fileWall->isopen());
  ofstream& fout = m_fileWall->get();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get current TRS
  SafePtr<TopologicalRegionSet> faceTrs = getCurrentTRS();

  // get SpectralFDMethodData
  SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  SafePtr<SpectralFDMethod> spectralFDMethod = spaceMethod.d_castTo<SpectralFDMethod>();
  cf_assert(spectralFDMethod.isNotNull());
  SafePtr<SpectralFDMethodData> spectralFDData = spectralFDMethod->getData();

  // get bndFacesStartIdxs from SpectralFDMethodData
  map< std::string , vector< vector< CFuint > > >&
    bndFacesStartIdxsPerTRS = spectralFDData->getBndFacesStartIdxs();
  vector< vector< CFuint > > bndFacesStartIdxs = bndFacesStartIdxsPerTRS[faceTrs->getName()];

  // number of face orientations (should be the same for all TRs)
  cf_assert(bndFacesStartIdxs.size() != 0);
  const CFuint nbOrients = bndFacesStartIdxs[0].size()-1;

  // number of TRs
  const CFuint nbTRs = faceTrs->getNbTRs();
  cf_assert(bndFacesStartIdxs.size() == nbTRs);

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoData = m_faceBuilder->getDataGE();
  geoData.cellsTRS = cellTrs;
  geoData.facesTRS = faceTrs;
  geoData.isBoundary = true;

  // SET THE BOUNDARY CONDITION STATE COMPUTER
  m_bndFaceTermComputer->setBcStateComputer(m_bcStateComputer);

  // GET FACE MAPPED COORDINATES OF OUTPUT POINTS
  vector< SpectralFDElementData* >& sdLocalData = spectralFDData->getSDLocalData();
  SafePtr< vector< RealVector > > outputPntsMappedCoord = sdLocalData[0]->getFaceOutputPntFaceMappedCoords();

  // loop over TRs
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    // loop over different orientations
    for (m_orient = 0; m_orient < nbOrients; ++m_orient)
    {
      // start and stop index of the faces with this orientation
      const CFuint startFaceIdx = bndFacesStartIdxs[iTR][m_orient  ];
      const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][m_orient+1];

      // SET THE ORIENTATION OF THE FACES
      m_bndFaceTermComputer->setFaceOrientation(m_orient);
      m_bndFaceTermComputer->setPointSet(*outputPntsMappedCoord);

      // loop over faces with this orientation
      for (CFuint faceID = startFaceIdx; faceID < stopFaceIdx; ++faceID)
      {
        // build the face GeometricEntity
        geoData.idx = faceID;
        m_face = m_faceBuilder->buildGE();

        // get the neighbouring cell
        m_intCell = m_face->getNeighborGeo(0);

        // get the states in the neighbouring cell
        m_cellStates = m_intCell->getStates();

        // if cell is parallel updatable, compute the output data
        if ((*m_cellStates)[0]->isParUpdatable())
        {
          // COMPUTE DATA IN BOUNDARY FACE TERM COMPUTER
          computeFaceData();

          // COMPUTE STATES IN OUTPUT POINTS
          const vector< RealVector > outStates = m_bndFaceTermComputer->reconstructGivenPntsStates(*m_cellStates);

          // LOOP OVER OUTPUT POINTS
          const CFuint nbrOutPnts = outStates.size();
          cf_assert(nbrOutPnts == outputPntsMappedCoord->size());
          for (CFuint iPnt = 0; iPnt < nbrOutPnts; ++iPnt)
          {
            // compute coordinates of output point
            const RealVector coord = m_face->computeCoordFromMappedCoord((*outputPntsMappedCoord)[iPnt]);

            // dereference current state
            const RealVector& state = outStates[iPnt];

            // pressure
            const CFreal rho = state[0];
            const CFreal invRho = 1./rho;
            CFreal rhoK2 = 0.0;
            for (CFuint iDim = 0; iDim < dim; ++iDim)
            {
              rhoK2 += state[iDim+1]*state[iDim+1];
            }
            rhoK2 *= 0.5*invRho;
            const CFreal p = gammaMinus1*(state[nbEqsM1] - rhoK2);

            // mach
            const CFreal a2 = gamma*p*invRho;
            const CFreal Mach = sqrt(2.0*invRho*rhoK2/a2);

            // temperature
            const CFreal T = a2/(gamma*R);

            // dimensional values
            const CFreal TDim   = T   * m_updateVarSet->getModel()->getTempRef();
            const CFreal pDim   = p   * m_updateVarSet->getModel()->getPressRef();
            const CFreal rhoDim = rho * (m_updateVarSet->getModel()->getReferencePhysicalData())[EulerTerm::RHO];

            // pressure coefficient
            const CFreal Cp = (pDim - m_pInf) / (0.5*m_rhoInf*m_uInf*m_uInf);

            // Output to File
            for (CFuint iDim = 0; iDim < dim; ++iDim)
            {
              fout << coord[iDim]  << " ";
            }
            fout << m_alphadeg   << " "
                 << pDim         << " "
                 << Mach         << " "
                 << Cp           << " "
                 << TDim         << " "
                 << rhoDim       << "\n";
          }
        }

        // release the face
        m_faceBuilder->releaseGE();
      }
    }
  }

  m_fileWall->close();
}

//////////////////////////////////////////////////////////////////////////////

void EulerConsComputeAeroSpectralFD::computeAero()
{
  CFAUTOTRACE;

  // get some data needed further
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint nbEqsM1 = nbEqs - 1;
  const CFreal gamma = m_updateVarSet->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  // initialize Cp and surface
  RealVector Cp(dim);
  Cp = 0.;
  CFreal surface = 0.0;

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cellTrs = MeshDataStack::getActive()->getTrs("InnerCells");

  // get current TRS
  SafePtr<TopologicalRegionSet> faceTrs = getCurrentTRS();

  // get SpectralFDMethodData
  SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  SafePtr<SpectralFDMethod> spectralFDMethod = spaceMethod.d_castTo<SpectralFDMethod>();
  cf_assert(spectralFDMethod.isNotNull());
  SafePtr<SpectralFDMethodData> spectralFDData = spectralFDMethod->getData();

  // get bndFacesStartIdxs from SpectralFDMethodData
  map< std::string , vector< vector< CFuint > > >&
    bndFacesStartIdxsPerTRS = spectralFDData->getBndFacesStartIdxs();
  vector< vector< CFuint > > bndFacesStartIdxs = bndFacesStartIdxsPerTRS[faceTrs->getName()];

  // number of face orientations (should be the same for all TRs)
  cf_assert(bndFacesStartIdxs.size() != 0);
  const CFuint nbOrients = bndFacesStartIdxs[0].size()-1;

  // number of TRs
  const CFuint nbTRs = faceTrs->getNbTRs();
  cf_assert(bndFacesStartIdxs.size() == nbTRs);

  // get the geodata of the face builder and set the TRSs
  FaceToCellGEBuilder::GeoData& geoData = m_faceBuilder->getDataGE();
  geoData.cellsTRS = cellTrs;
  geoData.facesTRS = faceTrs;
  geoData.isBoundary = true;

  // SET THE BOUNDARY CONDITION STATE COMPUTER
  m_bndFaceTermComputer->setBcStateComputer(m_bcStateComputer);

  // GET FACE MAPPED COORDINATES OF OUTPUT POINTS
  vector< SpectralFDElementData* >& sdLocalData = spectralFDData->getSDLocalData();
  SafePtr< vector< vector< RealVector > > > faceNodeMappedCoords = sdLocalData[0]->getFaceNodeCoords();

  // create integrator object (assume P2 geometrical polynomials) for integration over faces
  // and get quadrature point face mapped coordinates and wheights
  const CFuint solPolyOrder = sdLocalData[0]->getPolyOrder();
  CFuint integratorOrder = 0;
  if (dim == 2)
  {
    integratorOrder = solPolyOrder+1;
  }
  else if (dim == 3)
  {
    integratorOrder = solPolyOrder+3;
  }
  TensorProductGaussIntegrator tpIntegrator(static_cast<CFDim>(dim-1),static_cast<CFPolyOrder::Type>(integratorOrder));
  vector< RealVector > qPntsFaceMapCoords = tpIntegrator.getQuadPntsMappedCoords();
  vector< CFreal >     qPntsWheights  = tpIntegrator.getQuadPntsWheights();

  // loop over TRs
  for (CFuint iTR = 0; iTR < nbTRs; ++iTR)
  {
    // loop over different orientations
    for (m_orient = 0; m_orient < nbOrients; ++m_orient)
    {
      // start and stop index of the faces with this orientation
      const CFuint startFaceIdx = bndFacesStartIdxs[iTR][m_orient  ];
      const CFuint stopFaceIdx  = bndFacesStartIdxs[iTR][m_orient+1];

      // SET THE ORIENTATION OF THE FACES
      m_bndFaceTermComputer->setFaceOrientation(m_orient);
      m_bndFaceTermComputer->setPointSet(qPntsFaceMapCoords);

      // loop over faces with this orientation
      for (CFuint faceID = startFaceIdx; faceID < stopFaceIdx; ++faceID)
      {
        // build the face GeometricEntity
        geoData.idx = faceID;
        m_face = m_faceBuilder->buildGE();

        // get the neighbouring cell
        m_intCell = m_face->getNeighborGeo(0);

        // get the states in the neighbouring cell
        m_cellStates = m_intCell->getStates();

        // if cell is parallel updatable, compute the output data
        if ((*m_cellStates)[0]->isParUpdatable())
        {
          // COMPUTE DATA IN BOUNDARY FACE TERM COMPUTER
          computeFaceData();

          // GET FACE JACOBIAN VECTORS IN THE QUADRATURE POINTS
          const vector<RealVector>& jacobVecPntSet = *m_bndFaceTermComputer->getFaceJacobPntSet();

          // COMPUTE STATES IN QUADRATURE POINTS
          const vector< RealVector > pntStates = m_bndFaceTermComputer->reconstructGivenPntsStates(*m_cellStates);

          // COMPUTE INTEGRAL OF PRESSURE AND FACE SURFACE
          const CFuint nbrPnts = pntStates.size();
          RealVector cpFace(dim);
          cpFace = 0;
          CFreal surfFace = 0.0;
          cf_assert(nbrPnts == jacobVecPntSet.size());
          cf_assert(nbrPnts == qPntsWheights.size());
          for (CFuint iPnt = 0; iPnt < nbrPnts; ++iPnt)
          {
            // dereference current state
            const RealVector& state = pntStates[iPnt];

            // pressure
            const CFreal rho = state[0];
            const CFreal invRho = 1./rho;
            CFreal rhoK2 = 0.0;
            for (CFuint iDim = 0; iDim < dim; ++iDim)
            {
              rhoK2 += state[iDim+1]*state[iDim+1];
            }
            rhoK2 *= 0.5*invRho;
            const CFreal p = gammaMinus1*(state[nbEqsM1] - rhoK2);

            // add contribution to face integral
            cpFace += qPntsWheights[iPnt]*jacobVecPntSet[iPnt]*p;

            // add contribution to face
            surfFace += qPntsWheights[iPnt]*jacobVecPntSet[iPnt].norm2();
          }

          // add to global Cp and surface
          Cp += cpFace;
          surface += surfFace;
        }

        // release the face
        m_faceBuilder->releaseGE();
      }
    }
  }

  // adimensionalize Cp
  Cp *= m_updateVarSet->getModel()->getPressRef()/(0.5 * m_rhoInf * m_uInf * m_uInf * surface);

  // divide by reference length
  m_xForceCoef = Cp[XX];
  m_yForceCoef = Cp[YY];
  if (dim == 3)
  {
    m_zForceCoef = Cp[ZZ];
  }

  // project Cp (this should be changed for 3D)
  m_lift = (sin(m_alpharad)*Cp[XX] + cos(m_alpharad)*Cp[YY]);
  m_drag = (cos(m_alpharad)*Cp[XX] - sin(m_alpharad)*Cp[YY]);

  //Output to file the coefficients
  updateOutputFileAero();
}

//////////////////////////////////////////////////////////////////////////////

void EulerConsComputeAeroSpectralFD::setup()
{
  CFAUTOTRACE;

  // set lift and drag to zero
  m_lift = 0.;
  m_drag = 0.;

  // get Euler var set
  m_updateVarSet = getMethodData().getUpdateVarSet().d_castTo<EulerVarSet>();

  // get SpectralFDMethodData
  // assume that just one space method is available
  SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  SafePtr<SpectralFDMethod> spectralFDMethod = spaceMethod.d_castTo<SpectralFDMethod>();
  cf_assert(spectralFDMethod.isNotNull());
  SafePtr<SpectralFDMethodData> spectralFDData = spectralFDMethod->getData();

  // get face builder
  m_faceBuilder = spectralFDData->getFaceBuilder();

  // get boundary face term computer
  m_bndFaceTermComputer = spectralFDData->getBndFaceTermComputer();

  // get current TRS name
  SafePtr<TopologicalRegionSet> faceTrs = getCurrentTRS();
  const std::string currTRSName = faceTrs->getName();

  // get corresponding BCStateComputer
  SafePtr< vector< SafePtr< BCStateComputer > > >
      bcStateComputers = spectralFDData->getBCStateComputers();
  const CFuint nbrBCs = bcStateComputers->size();
  bool bcFound = false;
  for (CFuint iBC = 0; iBC < nbrBCs && !bcFound; ++iBC)
  {
    SafePtr< vector< std::string > > bcTRSNames = (*bcStateComputers)[iBC]->getTRSNames();
    const CFuint nbrTRSs = bcTRSNames->size();
    for (CFuint iTRS = 0; iTRS < nbrTRSs; ++iTRS)
    {
      if (currTRSName == (*bcTRSNames)[iTRS])
      {
        m_bcStateComputer = (*bcStateComputers)[iBC];
        bcFound = true;
        break;
      }
    }
  }

  if (!(m_uInf > 0. && m_rhoInf > 0. && m_pInf > 0.))
  {
    throw Common::BadValueException(FromHere(),"Values at infinity are not set!!!");
  }

  prepareOutputFileAero();
}

//////////////////////////////////////////////////////////////////////////////

void EulerConsComputeAeroSpectralFD::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void EulerConsComputeAeroSpectralFD::prepareOutputFileWall()
{
  CFAUTOTRACE;

  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  cf_assert (!m_fileWall->isopen());

  using namespace boost::filesystem;
  path file = Environment::DirPaths::getInstance().getResultsDir() / path ( m_nameOutputFileWall );
  file = PathAppender::getInstance().appendAllInfo( file );

  ofstream& fout = m_fileWall->open(file);

  fout << "TITLE = Wall_Values" << "\n";
  if (dim == 2)
  {
    fout << "VARIABLES = x y Alpha Pressure Mach Cp Temperature Density" << "\n";
  }
  else
  {
    fout << "VARIABLES = x y z Alpha Pressure Mach Cp Temperature Density" << "\n";
  }
}

//////////////////////////////////////////////////////////////////////////////

void EulerConsComputeAeroSpectralFD::prepareOutputFileAero()
{
  CFAUTOTRACE;

  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  using namespace boost::filesystem;
  path fpath = Environment::DirPaths::getInstance().getResultsDir() / path ( m_nameOutputFileAero );
  fpath = PathAppender::getInstance().appendParallel( fpath );

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& convergenceFile = fhandle->open(fpath);

  convergenceFile << "TITLE = Aerodynamic_Coeficients"  << "\n";
  if (dim == 2)
  {
    convergenceFile << "VARIABLES = Iter PhysTime Alpha LiftCoef DragCoef FxCoef FyCoef" << "\n";
  }
  else
  {
    convergenceFile << "VARIABLES = Iter PhysTime Alpha LiftCoef DragCoef FxCoef FyCoef FzCoef" << "\n";
  }

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void EulerConsComputeAeroSpectralFD::updateOutputFileAero()
{
  CFAUTOTRACE;

  using namespace boost::filesystem;
  path fpath = Environment::DirPaths::getInstance().getResultsDir() / path ( m_nameOutputFileAero );
  fpath = PathAppender::getInstance().appendParallel( fpath );

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& convergenceFile = fhandle->open(fpath, ios_base::app);

  Common::SafePtr<SubSystemStatus> subSysStatus =
    SubSystemStatusStack::getActive();

  convergenceFile
  << subSysStatus->getNbIter()        << " "
  << subSysStatus->getCurrentTime()   << " "
  << m_alphadeg                       << " "
  << m_lift                           << " "
  << m_drag                           << " "
  << m_xForceCoef                     << " "
  << m_yForceCoef                     << " ";
  if (PhysicalModelStack::getActive()->getDim() == 3)
  {
    convergenceFile << m_zForceCoef;
  }
  convergenceFile << "\n";

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace AeroCoeff

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

