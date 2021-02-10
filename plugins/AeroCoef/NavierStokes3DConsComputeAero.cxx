#include <numeric>

#include "Common/PE.hh"

#include "Common/CFMultiMap.hh"
#include "Common/SafePtr.hh"

#include "MathTools/RealMatrix.hh"
#include "MathTools/MathConsts.hh"

#include "Environment/DirPaths.hh"

#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/PathAppender.hh"
#include "Common/BadValueException.hh"
#include "Common/ParserException.hh"
#include "Framework/MeshData.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "Framework/PhysicalModel.hh"

#include "AeroCoef/AeroCoefFS.hh"
#include "AeroCoef/NavierStokes3DConsComputeAero.hh"

#include "FluctSplit/FluctuationSplit.hh"
#include "FluctSplit/InwardNormalsData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::FluctSplit;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace AeroCoef {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NavierStokes3DConsComputeAero,
                      DataProcessingData,
                      AeroCoefFSModule>
aNavierStokes3DConsComputeAeroProvider("NavierStokes3DConsComputeAero");

/////////////////////////////////////////////////////////////////////////////

void NavierStokes3DConsComputeAero::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("Alpha","Definition of the function defining the angle between the flow and the x-axis.");
   options.addConfigOption< CFuint >("SaveRateAero","Rate for saving the output file with aerodynamic coefficients.");
   options.addConfigOption< std::string >("OutputFileWall","Name of Output File to write the wall values.");
   options.addConfigOption< CFuint >("SaveRateWall","Save Output File containing the wall values every...iterations.");
   options.addConfigOption< std::string >("OutputFileAero","Name of Output File to write the results.");
   options.addConfigOption< CFreal >("UInf","Velocity at infinity.");
   options.addConfigOption< CFreal >("RhoInf","Density at infinity.");
   options.addConfigOption< CFreal >("PInf","Pressure at infinity.");
   options.addConfigOption< CFreal >("Aref","Reference area.");
   options.addConfigOption< CFreal >("Lref","Reference length.");
   options.addConfigOption< std::vector<CFreal> >("Xref","Reference point x to which the pressure momentum is computed");
   options.addConfigOption< CFreal >("TInf","Temperarure at infinity.");
   options.addConfigOption< bool >("AppendTime","Append time to file name.");
   options.addConfigOption< bool >("AppendIter","Append Iteration# to file name.");
}

//////////////////////////////////////////////////////////////////////////////

NavierStokes3DConsComputeAero::NavierStokes3DConsComputeAero( const std::string& name) :
  DataProcessingCom(name),
  m_sockets(),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_faceNeighCell("faceNeighCell"),
  m_updateVarSet(),
  m_vars("t"),
  m_eval(0.0,1),
  m_alphadeg(0.),
  m_alpharad(0.),
  m_diffVar(CFNULL),
  m_fsData(CFNULL),
  m_sumFMpres(0.,2*DIM_3D),
  m_sumFMfric(0.,2*DIM_3D),
  m_sumFMtotal(0.,2*DIM_3D),
  m_trsdata(0)
{
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
  setParameter("UInf",&m_uInf);

  m_rhoInf = 0.;
  setParameter("RhoInf",&m_rhoInf);

  m_pInf = 0.;
  setParameter("PInf",&m_pInf);

  m_tInf = 0.;
  setParameter("TInf",&m_tInf);

  m_Aref = 1.;
  setParameter("Aref",&m_Aref);

  m_Lref = 1.;
  setParameter("Lref",&m_Lref);

  m_Xref.resize(3);
  m_Xref[XX] = 0.0;
  m_Xref[YY] = 0.0;
  m_Xref[ZZ] = 0.0;
  setParameter("Xref",&m_Xref);

  m_appendIter = false;
  setParameter("AppendIter",&m_appendIter);

  m_appendTime = false;
  setParameter("AppendTime",&m_appendTime);

}

//////////////////////////////////////////////////////////////////////////////

NavierStokes3DConsComputeAero::~NavierStokes3DConsComputeAero()
{
	if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
NavierStokes3DConsComputeAero::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_faceNeighCell);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
NavierStokes3DConsComputeAero::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = m_sockets.getAllSourceSockets();
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DConsComputeAero::setup()
{
  SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  SafePtr<FluctuationSplit> fs = spaceMethod.d_castTo<FluctuationSplit>();

  m_fsData = fs->getData();
  if (m_fsData->getDiffusiveVar().isNotNull())
    if (m_fsData->getDiffusiveVar()->getName() != "Null" )
    {
      m_diffVar = m_fsData->getDiffusiveVar().d_castTo<NavierStokesVarSet>();
    }

  const RealVector& refValues = m_updateVarSet->getModel()->getReferencePhysicalData();
  if (refValues[EulerTerm::V] == 0.   ||
      refValues[EulerTerm::RHO] == 0. ||
      refValues[EulerTerm::P] == 0.   ||
     (refValues[EulerTerm::VX] == 0. && refValues[EulerTerm::VY] == 0. && refValues[EulerTerm::VZ] == 0.))
  {
    std::string msg = std::string("Some needed reference values not set.");
    throw BadValueException (FromHere(),msg);
  }

  m_trsdata.resize(getTrsList().size());
  for (CFuint i=0; i<getTrsList().size(); i++)
  {
    m_trsdata[i].m_isOverlapInitialized=false;
    m_trsdata[i].m_overlapFilter.resize(0);
    m_trsdata[i].m_nNodesInGeo=0;
    m_trsdata[i].m_inverseNumbering.resize(0);
    m_trsdata[i].m_faceLocalNumbering.resize(0);
    m_trsdata[i].m_nUpdatableFaces=0;
    m_trsdata[i].m_nNodes=0;
  }

  prepareOutputFileAero();
/*
  // DEBUG
  DataHandle<CFreal> cellextras = m_sockets.getSocketSource<CFreal>("cellextras")->getDataHandle();
  DataHandle<CFreal> stateextras = m_sockets.getSocketSource<CFreal>("stateextras")->getDataHandle();
  DataHandle<State*,GLOBAL> states  = socket_states.getDataHandle();
  stateextras.resize(states.size());
  for (CFuint i = 0; i < stateextras.size(); ++i)  stateextras[i]=1.;
  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");
  cellextras.resize(5*cells->getLocalNbGeoEnts());
  for (CFuint i = 0; i < cellextras.size(); ++i)  cellextras[i]=0.;
*/
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DConsComputeAero::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);

  // get the update varset
  m_updateVarSet = getMethodData().getUpdateVarSet().d_castTo<Physics::NavierStokes::EulerVarSet>();

  // Create the socket source for the boundaryNormals for each TRS to which this command applies
  // This is because we have different normals in FiniteVolume and FluctSplit / FiniteElement
  const CFuint nbTrs = getTrsNames().size();

  for (CFuint iTrs = 0; iTrs < nbTrs; ++iTrs)
  {
    std::string socketName = getTrsName(iTrs) + "-boundaryNormals";
    m_sockets.createSocketSource<const CFreal*>(socketName);
  }
/*
  // DEBUG
  m_sockets.createSocketSource<CFreal>("cellextras");
  m_sockets.createSocketSource<CFreal>("stateextras");
*/
  // set the function that computes the alpha
  try
  {
    setFunction();
  }
  catch (Common::Exception& e)
  {
    CFout << e.what() << "\n";
    throw;
  }

  // feedback from reference values
  CFLog(WARN, "PInf   : " << m_pInf << "\n");
  CFLog(WARN, "RhoInf : " << m_rhoInf << "\n");
  CFLog(WARN, "UInf   : " << m_uInf << "\n");
  CFLog(WARN, "Aref   : " << m_Aref << "\n");
  CFLog(WARN, "Lref   : " << m_Lref << "\n");
  CFLog(WARN, "Xref   : " << m_Xref[0] << " " << m_Xref[1] << " " << m_Xref[2] << "\n");

  if (MathChecks::isZero(m_pInf))   throw BadValueException (FromHere(),"pInf, the pressure at infinity, is zero");
  if (MathChecks::isZero(m_rhoInf)) throw BadValueException (FromHere(),"rhoInf, the density at infinity, is zero");
  if (MathChecks::isZero(m_uInf))   throw BadValueException (FromHere(),"uInf, the velocity at infinity, is zero");

}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DConsComputeAero::setFunction()
{
  // some sanity checks
  std::vector<std::string> m_functionDef = Common::StringOps::getWords(m_function);
  if(m_functionDef.size() != 1)
   throw BadValueException (FromHere(),"NavierStokes3DConsComputeAero::setFuntion(): alpha function wrongly defined.");

  m_functionParser.Parse(m_function, m_vars);

  if (m_functionParser.ErrorMsg() != 0)
  {
    std::string msg("ParseError in CFL::setFuntion(): ");
    msg += std::string(m_functionParser.ErrorMsg());
    msg += " Function: " + m_function;
    msg += " Vars: "     + m_vars;
    throw Common::ParserException (FromHere(),msg);
  }
}


//////////////////////////////////////////////////////////////////////////////

// the sort comparator
bool nodeComparator(CFuint* a, CFuint* b) {
  if (a[0]!=b[0]) return (a[0]<b[0]);
  for (int i=1; i<(const int)(a[0]+1); ++i)
    if (a[i]!=b[i]) return (a[i]<b[i]);
  return false;
}

void NavierStokes3DConsComputeAero::executeOnTrs()
{
  /// AL: this needs to be adapted to work with "long long int"
  const std::string nsp = getMethodData().getNamespace();
  
  CFAUTOTRACE;
  CFuint currtrsid=getCurrentTrsID();
  Common::SafePtr<Common::ConnectivityTable<CFuint> > g2n=getCurrentTRS()->getGeo2NodesConn();

  // if very first time, filling the trs-local renumbering of TRSDATA,
  //    CFuint m_nNodesInGeo;
  //    std::vector<CFint> m_faceLocalNumbering;
  //    std::vector<bool> m_useCoords
  // m_isOverlapInitialized will be taken care later
  // renumbering keeps the order of the nodes, so a linear for loop by testing with m_useCoords will dump coordinates accordingly
  if ((!m_trsdata[currtrsid].m_isOverlapInitialized)&&(g2n->nbRows()!=0))
  {

    // get connectivity, and find maxid for upper bound of renumber array
    if (g2n->size()!=g2n->nbCols(0)*g2n->nbRows()) throw Common::NotImplementedException(FromHere(),"Sorry, NavierStokes3DConsComputeAero only supports TRSs with same type elements.");
    m_trsdata[currtrsid].m_nNodesInGeo=g2n->nbCols(0); // !!!! checking first element' for nb nodes in an element, POTENTIALLY DANGEROUS !!!
    m_trsdata[currtrsid].m_faceLocalNumbering.reserve(g2n->size());
    CFuint maxID=0;

    for(CFuint i=0; i<g2n->nbRows(); ++i)
      for(CFuint j=0; j<g2n->nbCols(i); ++j)
        maxID=maxID>(*g2n)(i,j)?maxID:(*g2n)(i,j);
    m_trsdata[currtrsid].m_inverseNumbering.resize(maxID+1,-1);

    // compute mapping of face trs local numbering
    for(CFuint i=0; i<g2n->nbRows(); ++i)
      for(CFuint j=0; j<g2n->nbCols(i); ++j)
        if (m_trsdata[currtrsid].m_inverseNumbering[(*g2n)(i,j)]==-1)
          m_trsdata[currtrsid].m_inverseNumbering[(*g2n)(i,j)]=1;
    CFint currid=0;
    for(CFuint i=0; i<m_trsdata[currtrsid].m_inverseNumbering.size(); ++i)
      if (m_trsdata[currtrsid].m_inverseNumbering[i]!=-1)
        m_trsdata[currtrsid].m_inverseNumbering[i]=currid++;
    m_trsdata[currtrsid].m_nNodes=currid+1;

    // and apply
    for(CFuint i=0; i<g2n->nbRows(); ++i)
      for(CFuint j=0; j<g2n->nbCols(i); ++j)
        m_trsdata[currtrsid].m_faceLocalNumbering.push_back(m_trsdata[currtrsid].m_inverseNumbering[(*g2n)(i,j)]);
  }

  // the very first time initialize
  //    bool m_isOverlapInitialized;
  //    std::vector<bool> m_overlapFilter;
  if (PE::GetPE().GetProcessorCount(nsp)==1)
  {
    // fast exit on serial
    if (!m_trsdata[currtrsid].m_isOverlapInitialized)
    {
      m_trsdata[currtrsid].m_isOverlapInitialized=true;
      m_trsdata[currtrsid].m_overlapFilter.assign(getCurrentTRS()->getLocalNbGeoEnts(),true);
      m_trsdata[currtrsid].m_nUpdatableFaces=getCurrentTRS()->getLocalNbGeoEnts();
    }
  } else if (!m_trsdata[currtrsid].m_isOverlapInitialized)
  {

    // builder for standard TRS GeometricEntity's that will be used for the faces
    Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilderFace;
    geoBuilderFace.setup();
    StdTrsGeoBuilder::GeoData& geoDataFace = geoBuilderFace.getDataGE();
    geoDataFace.trs = getCurrentTRS();
    CFuint nbFaces = getCurrentTRS()->getLocalNbGeoEnts();
    const CFuint irank=PE::GetPE().GetRank(nsp);

    // the sorter data, one item is: nbnodes + node gids + process number + face local id
    std::vector<CFuint> rproc(0);
    rproc.reserve(nbFaces*13+1); // avoid memory issues on mpi communication
    std::vector<CFuint> rglob(0);
    rglob.reserve(1); // avoid memory issues on mpi communication

    // loop on faces to fill
    for (CFuint iFace = 0; iFace < nbFaces; ++iFace)
    {
      geoDataFace.idx = iFace;
      GeometricEntity & currFace = *geoBuilderFace.buildGE();
      const CFuint nbGeoNodes=currFace.nbNodes();
      rproc.push_back(nbGeoNodes);
      for(CFuint iNode=0 ; iNode<nbGeoNodes; ++iNode)
        rproc.push_back(currFace.getNode(iNode)->getGlobalID());
      rproc.push_back(irank);
      rproc.push_back(iFace);
      geoBuilderFace.releaseGE();
      sort(rproc.end()-(nbGeoNodes+2),rproc.end()-2);
    }

    // collect, prepare and sort
    std::vector<int> nb(PE::GetPE().GetProcessorCount(nsp),0);
    int nbproc=(int)rproc.size();
    MPI_Gather(&nbproc, 1, MPI_INT, &nb[0], 1, MPI_INT, 0, PE::GetPE().GetCommunicator(nsp));
    int totnb=std::accumulate(nb.begin(),nb.end(),0);
    if (irank==0) rglob.resize(totnb);
    std::vector<int> displs(PE::GetPE().GetProcessorCount(nsp),0);
    for(int i=1; i<(const int)displs.size(); i++) displs[i]=displs[i-1]+nb[i-1];
    MPI_Gatherv(&rproc[0],nbproc,MPI_INT,&rglob[0],&nb[0], &displs[0],MPI_INT,0, PE::GetPE().GetCommunicator(nsp));
    std::vector<CFuint*> pglob(0);
    int npglob=0;
    for(int i=0; i<totnb; i+=rglob[i]+3) ++npglob;
    pglob.reserve(npglob+1);  // avoid memory issues on mpi communication
    for(int i=0; i<totnb; i+=rglob[i]+3) pglob.push_back(&rglob[i]);
    std::sort(pglob.begin(),pglob.end(),nodeComparator);

    // build list of duplicates and sort by rank (bruteforce, assuming small)
    std::vector<std::pair<CFuint,CFuint> > dup(0);
    for(int i=1; i<(const int)(pglob.size()); ++i)
    {
      if (pglob[i-1][0]!=pglob[i][0]) continue;
      CFuint j=1;
      for(; j<(const CFuint)(pglob[i][0]+1); ++j)
        if (pglob[i-1][j]!=pglob[i][j])
          break;
      if (j!=pglob[i][0]+1)
        continue;
      dup.push_back(std::pair<CFuint,CFuint>(pglob[i][pglob[i][0]+1],pglob[i][pglob[i][0]+2]));
    }
    std::sort(dup.begin(),dup.end());
    std::vector<int> duplid(0);
    duplid.reserve(dup.size());
    for(int i=0; i<(const int)(dup.size()); ++i) duplid.push_back(dup[i].second);

    // distribute back from process 0
    nb.assign(PE::GetPE().GetProcessorCount(nsp),0);
    for (int i=0; i<(const int)dup.size(); i++) nb[dup[i].first]++;
    displs.assign(PE::GetPE().GetProcessorCount(nsp),0);
    for(int i=1; i<(const int)displs.size(); i++) displs[i]=displs[i-1]+nb[i-1];
    MPI_Scatter(&nb[0], 1, MPI_INT, &nbproc, 1, MPI_INT, 0, PE::GetPE().GetCommunicator(nsp));
    std::vector<int> localdup(nbproc,0);
    MPI_Scatterv(&duplid[0],&nb[0],&displs[0],MPI_INT,&localdup[0],nbproc,MPI_INT,0, PE::GetPE().GetCommunicator(nsp));

    // finally mark duplicates as false
    m_trsdata[currtrsid].m_overlapFilter.assign(nbFaces,true);
    for (int i=0; i<(const int)localdup.size(); ++i) m_trsdata[currtrsid].m_overlapFilter[localdup[i]]=false;
    m_trsdata[currtrsid].m_isOverlapInitialized=true;

    // test of active number of faces
    int nLocalUpdatableFace=0;
    for (int i=0; i<(const int)m_trsdata[currtrsid].m_overlapFilter.size(); ++i)
      if (m_trsdata[currtrsid].m_overlapFilter[i])
        ++nLocalUpdatableFace;
    int nUpdatableFace=0;
    MPI_Allreduce(&nLocalUpdatableFace,  &nUpdatableFace,  1, MPI_INT, MPI_SUM, PE::GetPE().GetCommunicator(nsp));
    CFLog(WARN,"Number of wall TRS's faces after filtering: " << nUpdatableFace << "\n" );
    m_trsdata[currtrsid].m_nUpdatableFaces=nLocalUpdatableFace;
  }

  // compute coefficients if necessary
  Common::SafePtr<SubSystemStatus> ssys_status = SubSystemStatusStack::getActive();
  const CFuint iter = ssys_status->getNbIter();
  if((!(iter % m_saveRateAero))||(!(iter % m_saveRateWall))) { computeWall(!(iter % m_saveRateWall),!(iter % m_saveRateAero)); }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DConsComputeAero::computeCellVelocityAndFacePressureAndDynamicViscosity( COOLFluiD::Framework::GeometricEntity& currFace, COOLFluiD::Framework::GeometricEntity& neighborCell, std::vector<RealVector>& vel, RealVector& pres, CFreal& mu )
{
  // setup
  vector<Framework::State*>& cstates= (*(neighborCell.getStates()));
  vector<Framework::State*>& fstates= (*(currFace.getStates()));
  vel.resize(cstates.size(),RealVector(0.,DIM_3D));
  pres.resize(fstates.size());

  // fill velocity
  for(CFuint iState=0; iState<cstates.size(); ++iState)
  {
    const CFreal invRho=1./(*cstates[iState])[0];
    vel[iState][0]=(*cstates[iState])[1]*invRho;
    vel[iState][1]=(*cstates[iState])[2]*invRho;
    vel[iState][2]=(*cstates[iState])[3]*invRho;
  }

  // fill pressure
  const CFreal gamma = m_updateVarSet->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;
  for(CFuint iState=0; iState<fstates.size(); ++iState)
  {
    const CFreal rho  = (*cstates[iState])[0];
    const CFreal rhoU = (*cstates[iState])[1];
    const CFreal rhoV = (*cstates[iState])[2];
    const CFreal rhoW = (*cstates[iState])[3];
    const CFreal rhoE = (*cstates[iState])[4];
    const CFreal halfRhoVel2 = 0.5/rho*(rhoU*rhoU+rhoV*rhoV+rhoW*rhoW);
    pres[iState] = gammaMinus1*(rhoE - halfRhoVel2);
  }

  // set dynamic viscosity
  mu = 0.;
  if (m_diffVar.isNotNull()) mu = m_diffVar->getCurrDynViscosity();
}

//////////////////////////////////////////////////////////////////////////////

std::vector< RealVector > NavierStokes3DConsComputeAero::computeWallStateUnitNormals(GeometricEntity & currFace, CFuint iFace, DataHandle<const CFreal*> boundaryNormals)
{
  std::vector< RealVector > retvec(0);
  for(CFuint iState = 0; iState < currFace.nbStates(); iState++){
    State* istate=(*currFace.getStates())[iState];
    RealVector istatenorm=currFace.computeCellNormal(*(istate->getCoordinates().getData()));
    istatenorm.normalize();
    retvec.push_back(istatenorm);
  }
  return retvec;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< RealMatrix > NavierStokes3DConsComputeAero::computeWallStateVelocityGradients( GeometricEntity & currFace, GeometricEntity& neighborCell, std::vector<RealVector>& vel )
{
  // return data
  std::vector< RealMatrix > retvec(0);

  // get face/cell states and face coordinates
  vector<Framework::State*>& cstates= (*(neighborCell.getStates()));
  const CFuint nbStatesInCell=cstates.size();
  vector<Framework::State*>& fstates= (*(currFace.getStates()));
  const CFuint nbStatesInFace=fstates.size();
  vector<RealVector> fcoords(0);
  for(CFuint iState=0; iState<nbStatesInFace; iState++) fcoords.push_back(*(fstates[iState]->getCoordinates().getData()));
  //cout << "FCOORDS(" << fcoords.size() << ") " << fcoords[0] << " " << fcoords[1] << " " << fcoords[2] << "\n" << flush;

  // Compute the shape function gradient at the projected point
  const std::vector<RealMatrix> cellGradients= neighborCell.computeSolutionShapeFunctionGradients(fcoords);

  // compute gradients
  for (CFuint iState = 0; iState < nbStatesInFace; ++iState){
    retvec.push_back(RealMatrix(3,3));
    retvec[iState]=0.;
    for (CFuint iNode = 0; iNode < nbStatesInCell ; ++iNode)
      for (CFuint iDim = 0; iDim < DIM_3D ; ++iDim)
        for (CFuint jDim = 0; jDim < DIM_3D ; ++jDim)
          (retvec[iState])(jDim,iDim) += (cellGradients[iState])(iNode,iDim) *vel[iNode][jDim];
  }

  return retvec;
}

//////////////////////////////////////////////////////////////////////////////

RealVector NavierStokes3DConsComputeAero::computeForceAndMoment(COOLFluiD::Framework::GeometricEntity& currFace, std::vector<RealVector> tau)
{
  // from stress to forces & moments, involving integration on the wall face
  // first 3 are force, second 3 are moment
  // this is real force and moment, it is non-dimensionalized with free stream values somewhere later
  // in this form it only works for P1P1 and P1P2 tetrahedrons
  RealVector fm(0.,2*DIM_3D);
  if (currFace.getSolutionShapeFunctionOrder()==CFPolyOrder::ORDER1)
  {
    fm[0]=(tau[0][0]+tau[1][0]+tau[2][0]);
    fm[1]=(tau[0][1]+tau[1][1]+tau[2][1]);
    fm[2]=(tau[0][2]+tau[1][2]+tau[2][2]);
    for (CFuint iState=0; iState<currFace.nbStates(); iState++)
    {
      CFreal r0=(*(currFace.getNode(iState)))[0]-m_Xref[0];
      CFreal r1=(*(currFace.getNode(iState)))[1]-m_Xref[1];
      CFreal r2=(*(currFace.getNode(iState)))[2]-m_Xref[2];
      fm[3]+=tau[iState][2]*r1-tau[iState][1]*r2;
      fm[4]+=tau[iState][0]*r2-tau[iState][2]*r0;
      fm[5]+=tau[iState][1]*r0-tau[iState][0]*r1;
    }
    const CFreal dx01=(*(currFace.getNode(1)))[0]-(*(currFace.getNode(0)))[0];
    const CFreal dy01=(*(currFace.getNode(1)))[1]-(*(currFace.getNode(0)))[1];
    const CFreal dz01=(*(currFace.getNode(1)))[2]-(*(currFace.getNode(0)))[2];
    const CFreal dx02=(*(currFace.getNode(2)))[0]-(*(currFace.getNode(0)))[0];
    const CFreal dy02=(*(currFace.getNode(2)))[1]-(*(currFace.getNode(0)))[1];
    const CFreal dz02=(*(currFace.getNode(2)))[2]-(*(currFace.getNode(0)))[2];
    const CFreal area=0.5*fabs(+dy01*dz02-dy02*dz01
                               -dx01*dz02+dx02*dz01
                               +dx01*dy02-dx02*dy01);
    fm*=area/3.;
  }
  if (currFace.getSolutionShapeFunctionOrder()==CFPolyOrder::ORDER2)
  {
    fm[0]=(tau[3][0]+tau[4][0]+tau[5][0]);
    fm[1]=(tau[3][1]+tau[4][1]+tau[5][1]);
    fm[2]=(tau[3][2]+tau[4][2]+tau[5][2]);
    for (CFuint iState=3; iState<currFace.nbStates(); iState++)
    {
      CFreal r0=(*(currFace.getNode(iState)))[0]-m_Xref[0];
      CFreal r1=(*(currFace.getNode(iState)))[1]-m_Xref[1];
      CFreal r2=(*(currFace.getNode(iState)))[2]-m_Xref[2];
      fm[3]+=tau[iState][2]*r1-tau[iState][1]*r2;
      fm[4]+=tau[iState][0]*r2-tau[iState][2]*r0;
      fm[5]+=tau[iState][1]*r0-tau[iState][0]*r1;
    }
    const CFreal dx01=(*(currFace.getNode(1)))[0]-(*(currFace.getNode(0)))[0];
    const CFreal dy01=(*(currFace.getNode(1)))[1]-(*(currFace.getNode(0)))[1];
    const CFreal dz01=(*(currFace.getNode(1)))[2]-(*(currFace.getNode(0)))[2];
    const CFreal dx02=(*(currFace.getNode(2)))[0]-(*(currFace.getNode(0)))[0];
    const CFreal dy02=(*(currFace.getNode(2)))[1]-(*(currFace.getNode(0)))[1];
    const CFreal dz02=(*(currFace.getNode(2)))[2]-(*(currFace.getNode(0)))[2];
    const CFreal area=0.5*fabs(+dy01*dz02-dy02*dz01
                               -dx01*dz02+dx02*dz01
                               +dx01*dy02-dx02*dy01);
    fm*=area/3.;
  }
  return fm;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DConsComputeAero::computeWall( bool saveWall, bool saveAero)
{
  CFAUTOTRACE;
  
  const std::string nsp = getMethodData().getNamespace();
  
  const CFuint irank=PE::GetPE().GetRank(nsp);
  const CFuint nproc=PE::GetPE().GetProcessorCount(nsp);
  const bool writeWall=(saveWall)&&(getCurrentTRS()->getLocalNbGeoEnts()!=0);
  const bool writeAero=(saveAero)&&(irank==0);
  CFuint currtrsid=getCurrentTrsID();

  // prepare wall output file
  if (writeWall) prepareOutputFileWall();

  // get the trs and datahandles
  const std::string TRSName = getCurrentTRS()->getName();
  const std::string socketName = TRSName + "-boundaryNormals";
  DataHandle<const CFreal*> boundaryNormals = m_sockets.getSocketSource<const CFreal*>(socketName)->getDataHandle();
  DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > > faceNeighCell = socket_faceNeighCell.getDataHandle();
/*
  // DEBUG
  DataHandle<CFreal> cellextras = m_sockets.getSocketSource<CFreal>("cellextras")->getDataHandle();
  DataHandle<CFreal> stateextras = m_sockets.getSocketSource<CFreal>("stateextras")->getDataHandle();
*/
  // builder for standard TRS GeometricEntity's that will be used for the faces
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilderFace;
  geoBuilderFace.setup();
  StdTrsGeoBuilder::GeoData& geoDataFace = geoBuilderFace.getDataGE();
  geoDataFace.trs = getCurrentTRS();
  const CFuint nbFaces = getCurrentTRS()->getLocalNbGeoEnts();
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilderCell;
  geoBuilderCell.setup();
  StdTrsGeoBuilder::GeoData& geoDataCell = geoBuilderCell.getDataGE();

  // data for the wall surface file
  std::vector<std::vector<CFreal> > wall_values(29, std::vector<CFreal>(0) );
  for (CFuint i=0; i<wall_values.size(); ++i) wall_values[i].reserve(m_trsdata[currtrsid].m_nUpdatableFaces);
  std::vector<CFreal> wall_coords(m_trsdata[currtrsid].m_nNodes*DIM_3D);

  // this is the sum
  RealVector sumFMfric(0.,2*DIM_3D);
  RealVector sumFMpres(0.,2*DIM_3D);
  RealVector sumFMtotal(0.,2*DIM_3D);
    
  for(CFuint iProc=0; iProc<nproc; ++iProc){
    MPI_Barrier( PE::GetPE().GetCommunicator(nsp));
    if(iProc==irank){

      // loop on the faces
      for (CFuint iFace = 0; iFace < nbFaces; ++iFace)
      {
        // build boundary face geoent
        geoDataFace.idx = iFace;
        GeometricEntity & currFace = *geoBuilderFace.buildGE();

        // normal at the states
        std::vector<RealVector> normal=computeWallStateUnitNormals(currFace,iFace,boundaryNormals);
        //cout << "NORMAL(" << normal.size() << ") " << normal[0] << "   " << normal[1] << "   " << normal[2] << "\n" << flush;

        // bnd face centroid coordinates
        RealVector centroid(0.,DIM_3D);
        for(int iState=0; iState<(const int)currFace.nbStates(); ++iState)
          for(int iDim=0; iDim<DIM_3D; ++iDim)
            centroid[iDim]+=currFace.getState(iState)->getCoordinates()[iDim];
        centroid/=(CFreal)currFace.nbStates();

        // compute gradient at walls (at bnd face states)
        geoDataCell.idx = faceNeighCell[iFace].first;
        geoDataCell.trs = faceNeighCell[iFace].third;
        GeometricEntity& neighborCell = *geoBuilderCell.buildGE();
/*
        // DEBUG
        cellextras[5*neighborCell.getID()+0]=1.;
        cellextras[5*neighborCell.getID()+1]=normal[0][0];
        cellextras[5*neighborCell.getID()+2]=normal[1][1];
        cellextras[5*neighborCell.getID()+3]=normal[2][2];
        const CFreal dx01=(*(currFace.getNode(1)))[0]-(*(currFace.getNode(0)))[0];
        const CFreal dy01=(*(currFace.getNode(1)))[1]-(*(currFace.getNode(0)))[1];
        const CFreal dz01=(*(currFace.getNode(1)))[2]-(*(currFace.getNode(0)))[2];
        const CFreal dx02=(*(currFace.getNode(2)))[0]-(*(currFace.getNode(0)))[0];
        const CFreal dy02=(*(currFace.getNode(2)))[1]-(*(currFace.getNode(0)))[1];
        const CFreal dz02=(*(currFace.getNode(2)))[2]-(*(currFace.getNode(0)))[2];
        cellextras[5*neighborCell.getID()+4]=0.5*fabs(+dy01*dz02-dy02*dz01
                                                      -dx01*dz02+dx02*dz01
                                                      +dx01*dy02-dx02*dy01);
        for(int iState=0; iState<currFace.nbStates(); ++iState)
          stateextras[currFace.getState(iState)->getLocalID()]=0.;
*/
        // compute variables
        std::vector<RealVector> vel;
        RealVector pres(0.,DIM_3D);
        CFreal mu;
        computeCellVelocityAndFacePressureAndDynamicViscosity(currFace,neighborCell,vel,pres,mu);
        pres-=m_pInf;
        //cout << "VEL(" << vel.size() << ") " << vel[0] << "   " << vel[1] << "   " << vel[2] << "   " << vel[3] << "\n" << flush;
        //cout << "PRES(" << pres.size() << ") " << pres << "\n" << flush;
        //cout << "MU " << mu << "\n" << flush;
        std::vector<RealMatrix> velgrad=computeWallStateVelocityGradients(currFace,neighborCell,vel);
        //cout << "VELGRAD(" << velgrad.size() << ") " << velgrad[0] << "   " << velgrad[1] << "   " << velgrad[2] << "\n" << flush;
        std::vector<RealVector> taufric(0);
        std::vector<RealVector> taupres(0);
        for(CFuint iState=0; iState<currFace.nbStates(); ++iState)
        {
          RealMatrix stress(DIM_3D,DIM_3D);
          velgrad[iState].transpose(stress);
          stress+=velgrad[iState];
          stress*=mu;
          taufric.push_back(stress*normal[iState]);
          CFreal permag=0.;
          for(int idim=0; idim<DIM_3D; idim++)
            permag+=taufric[iState][idim]*normal[iState][idim];
          taufric[iState]-=permag*normal[iState];
          taupres.push_back(pres[iState]*normal[iState]);
          taupres[iState]*=-1;
        }

        // integrate the friction/pressure and return forces and moments
        RealVector FMfric=computeForceAndMoment(currFace,taufric);
        RealVector FMpres=computeForceAndMoment(currFace,taupres);
        RealVector FMtotal=FMfric+FMpres;
        if (m_trsdata[currtrsid].m_overlapFilter[iFace])
        {
          sumFMfric+= FMfric;
          sumFMpres+= FMpres;
          sumFMtotal+=FMtotal;
        }

        // collect surface data
        if ((m_trsdata[currtrsid].m_overlapFilter[iFace])&&(writeWall)){

          // data
          int ient=0;
          wall_values[ient++].push_back(normal[0][0]);
          wall_values[ient++].push_back(normal[0][1]);
          wall_values[ient++].push_back(normal[0][2]);
          wall_values[ient++].push_back(FMtotal[0]);
          wall_values[ient++].push_back(FMtotal[1]);
          wall_values[ient++].push_back(FMtotal[2]);
          wall_values[ient++].push_back(FMtotal[3]);
          wall_values[ient++].push_back(FMtotal[4]);
          wall_values[ient++].push_back(FMtotal[5]);
          wall_values[ient++].push_back(FMfric[0]);
          wall_values[ient++].push_back(FMfric[1]);
          wall_values[ient++].push_back(FMfric[2]);
          wall_values[ient++].push_back(FMfric[3]);
          wall_values[ient++].push_back(FMfric[4]);
          wall_values[ient++].push_back(FMfric[5]);
          wall_values[ient++].push_back(FMpres[0]);
          wall_values[ient++].push_back(FMpres[1]);
          wall_values[ient++].push_back(FMpres[2]);
          wall_values[ient++].push_back(FMpres[3]);
          wall_values[ient++].push_back(FMpres[4]);
          wall_values[ient++].push_back(FMpres[5]);
          wall_values[ient++].push_back(m_rhoInf);
          wall_values[ient++].push_back(m_uInf);
          wall_values[ient++].push_back(m_pInf);
          wall_values[ient++].push_back(m_Aref);
          wall_values[ient++].push_back(m_Lref);
          wall_values[ient++].push_back(m_Xref[0]);
          wall_values[ient++].push_back(m_Xref[1]);
          wall_values[ient++].push_back(m_Xref[2]);

          // coords
          for(int iNode=0; iNode<(const int)currFace.nbNodes(); ++iNode)
            if (m_trsdata[currtrsid].m_inverseNumbering[currFace.getNode(iNode)->getLocalID()]!=-1)
              for(int iDim=0; iDim<DIM_3D; ++iDim)
                wall_coords[m_trsdata[currtrsid].m_inverseNumbering[currFace.getNode(iNode)->getLocalID()]*3+iDim]=(*(currFace.getNode(iNode)))[iDim];
        }
/*
        // write surface file
        if ((m_trsdata[currtrsid].m_overlapFilter[iFace])&&(writeWall)){
          (*fout) << centroid  << " "
                  << normal[0] << " "
                  << FMtotal   << " "
                  << FMfric    << " "
                  << FMpres    << " "
                  << m_rhoInf  << " "
                  << m_uInf    << " "
                  << m_pInf    << " "
                  << m_Aref    << " "
                  << m_Lref    << " "
                  << m_Xref[0] << " " << m_Xref[1] << " " << m_Xref[2] << " "
                  << "\n" << flush;
        }
*/
        // get rid of stuff
        geoBuilderFace.releaseGE();
        geoBuilderCell.releaseGE();

      }

      if (writeWall){

        // opening wall file on process
        std::string walltrsname(getCurrentTRS()->getName()+"-"+m_nameOutputFileWall);
        boost::filesystem::path fpath = Environment::DirPaths::getInstance().getResultsDir() / boost::filesystem::path (walltrsname);
        fpath = PathAppender::getInstance().appendAllInfo(fpath, m_appendIter, m_appendTime, true );
        SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
        ofstream* fout = &(fhandle->open(fpath));

        // print header
        (*fout) << "TITLE  =  \"Values on the Wall\"" << "\n";
        (*fout) << "VARIABLES = x y z nx ny nz Fx Fy Fz Mx My Mz Ffx Ffy Ffz Mfx Mfy Mfz Fpx Fpy Fpz Mpx Mpy Mpz Rhoinf Uinf Pinf Aref Lref Lmxref Lmyref Lmzref" << "\n" << flush;
        (*fout) << "ZONE T=\"" << getCurrentTRS()->getName() << "\", DATAPACKING=BLOCK, N=" << m_trsdata[currtrsid].m_nNodes << ", E=" << m_trsdata[currtrsid].m_nUpdatableFaces << ", ZONETYPE=FETRIANGLE, VARLOCATION=([4-32]=CELLCENTERED)\n" << flush;

        // write coordinates
        for (int idim=0; idim<DIM_3D; ++idim)
        {
          for (CFuint i=0; i<m_trsdata[currtrsid].m_nNodes; ++i)
            { (*fout) << " " << wall_coords[i*DIM_3D+idim]; if ((i+1)%10==0) (*fout) << "\n"; }
          if (m_trsdata[currtrsid].m_nNodes%10!=0) (*fout) << "\n";
        }

        // write data
        for (CFuint idat=0; idat<wall_values.size(); ++idat)
        {
          int i=0;
          for (CFuint ifil=0; ifil<m_trsdata[currtrsid].m_overlapFilter.size(); ++ifil)
            if (m_trsdata[currtrsid].m_overlapFilter[ifil])
              { (*fout) << " " << wall_values[idat][i]; if ((i+1)%10==0) (*fout) << "\n"; ++i; }
          if (m_trsdata[currtrsid].m_nUpdatableFaces%10!=0) (*fout) << "\n";
        }

        // write connectivity
        for (CFuint ifac=0; ifac<getCurrentTRS()->getLocalNbGeoEnts(); ++ifac)
          if (m_trsdata[currtrsid].m_overlapFilter[ifac])
          {
            for (CFuint i=0; i<m_trsdata[currtrsid].m_nNodesInGeo; ++i)
              (*fout) << " " << m_trsdata[currtrsid].m_faceLocalNumbering[ifac*m_trsdata[currtrsid].m_nNodesInGeo+i]+1;
            (*fout) << "\n";
          }

        fhandle->close();
      }

      MPI_Barrier( PE::GetPE().GetCommunicator(nsp));
    }
  }


  // compute sum and convert to coefficients
  MPI_Allreduce(&sumFMfric[0],  &m_sumFMfric[0],  2*DIM_3D, MPI_DOUBLE, MPI_SUM, PE::GetPE().GetCommunicator(nsp));
  MPI_Allreduce(&sumFMpres[0],  &m_sumFMpres[0],  2*DIM_3D, MPI_DOUBLE, MPI_SUM, PE::GetPE().GetCommunicator(nsp));
  MPI_Allreduce(&sumFMtotal[0], &m_sumFMtotal[0], 2*DIM_3D, MPI_DOUBLE, MPI_SUM, PE::GetPE().GetCommunicator(nsp));
  m_sumFMpres /=0.5*m_rhoInf*m_uInf*m_uInf*m_Aref;
  m_sumFMpres[3] /=m_Lref;
  m_sumFMpres[4] /=m_Lref;
  m_sumFMpres[5] /=m_Lref;
  m_sumFMfric /=0.5*m_rhoInf*m_uInf*m_uInf*m_Aref;
  m_sumFMfric[3] /=m_Lref;
  m_sumFMfric[4] /=m_Lref;
  m_sumFMfric[5] /=m_Lref;
  m_sumFMtotal/=0.5*m_rhoInf*m_uInf*m_uInf*m_Aref;
  m_sumFMtotal[3]/=m_Lref;
  m_sumFMtotal[4]/=m_Lref;
  m_sumFMtotal[5]/=m_Lref;

  if (writeAero) {
    cout << "Global aerodynamic coefficients: " << m_sumFMtotal << "    " << m_sumFMfric << "    " << m_sumFMpres << "\n" << flush;
    using namespace boost::filesystem;
    path fpath = Environment::DirPaths::getInstance().getResultsDir() / path ( m_nameOutputFileAero );
    fpath = PathAppender::getInstance().appendAllInfo(fpath, false, false, false);
    SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
    ofstream& convergenceFile = fhandle->open(fpath,ios_base::app);
    Common::SafePtr<SubSystemStatus> subSysStatus = SubSystemStatusStack::getActive();
    convergenceFile << subSysStatus->getNbIter()       << " "
                    << subSysStatus->getCurrentTime()  << " "
                    << m_sumFMtotal                    << " "
                    << m_sumFMfric                     << " "
                    << m_sumFMpres                     << " "
                    << m_rhoInf                        << " "
                    << m_uInf                          << " "
                    << m_pInf                          << " "
                    << m_Aref                          << " "
                    << m_Lref                          << " "
                    << m_Xref[0] << " " << m_Xref[1] << " " << m_Xref[2] << " "
                    << "\n" << flush;
    fhandle->close();
  }


}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DConsComputeAero::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DConsComputeAero::prepareOutputFileWall()
{
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokes3DConsComputeAero::prepareOutputFileAero()
{
  using namespace boost::filesystem;
  path fpath = Environment::DirPaths::getInstance().getResultsDir() / path ( m_nameOutputFileAero );
  fpath = PathAppender::getInstance().appendAllInfo(fpath, false, false, false);
  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& convergenceFile = fhandle->open(fpath);
  convergenceFile << "TITLE  =  \"Aerodynamic Coefficients\""  << "\n";
  convergenceFile << "VARIABLES = Iter PhysTime Cx Cy Cz Cmx Cmy Cmz Cfx Cfy Cfz Cmfx Cmfy Cmfz Cpx Cpy Cpz Cmpx Cmpy Cmpz Rhoinf Uinf Pinf Aref Lref Lmxref Lmyref Lmzref" << "\n" << flush;
  convergenceFile.close();
  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

