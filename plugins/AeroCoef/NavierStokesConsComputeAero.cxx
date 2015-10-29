#include "Common/PE.hh"
#include "Common/CFMultiMap.hh"
#include "Common/SafePtr.hh"
#include "MathTools/RealMatrix.hh"
#include "Environment/DirPaths.hh"

#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/PathAppender.hh"
#include "Common/BadValueException.hh"
#include "Common/ParserException.hh"
#include "Framework/GeometricEntity.hh"
#include "Framework/MeshData.hh"
#include "Framework/CellTrsGeoBuilder.hh"
#include "Framework/MapGeoToTrsAndIdx.hh"
#include "Framework/PhysicalModel.hh"

#include "AeroCoef/AeroCoefFS.hh"
#include "AeroCoef/NavierStokesConsComputeAero.hh"

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

MethodCommandProvider<NavierStokesConsComputeAero,
                      DataProcessingData,
                      AeroCoefFSModule>
aNavierStokesConsComputeAeroProvider("NavierStokesConsComputeAero");

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAero::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< std::string >("Alpha","Function defining the angle between the flow and the x-axis.");
   options.addConfigOption< std::string >("Beta","Function defining the angle between the flow and the y-axis.");
   options.addConfigOption< CFuint >("SaveRateAero","Rate for saving the output file with aerodynamic coeficients.");
   options.addConfigOption< std::string >("OutputFileWall","Name of Output File to write the wall values.");
   options.addConfigOption< CFuint >("SaveRateWall","Save Output File containing the wall values every...iterations.");
   options.addConfigOption< std::string >("OutputFileAero","Name of Output File to write the results.");
   options.addConfigOption< CFreal >("UInf","Velocity at infinity.");
   options.addConfigOption< CFreal >("RhoInf","Density at infinity.");
   options.addConfigOption< CFreal >("PInf","Pressure at infinity.");
   options.addConfigOption< std::vector<CFreal> >("Xref","Reference point x to which the pressure momentum is computed");
   options.addConfigOption< CFreal >("TInf","Temperarure at infinity.");
   options.addConfigOption< bool >("AppendTime","Append time to file name.");
   options.addConfigOption< bool >("AppendIter","Append Iteration# to file name.");

}

//////////////////////////////////////////////////////////////////////////////

NavierStokesConsComputeAero::NavierStokesConsComputeAero(const std::string& name) :
  DataProcessingCom(name),
  m_sockets(),
  socket_nodes("nodes"),
  socket_states("states"),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  m_updateVarSet(),
  m_cells(CFNULL),
  m_fsData(CFNULL),
  m_avValues(),
  m_diffVar(),
  m_vars("t"),
  m_eval(0.0,1),
  m_alphadeg(0.),
  m_betadeg(0.),
  m_dimState()
{
  m_fileWall = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();

  addConfigOptionsTo(this);

  m_nameOutputFileWall = "Wall";
  setParameter("OutputFileWall",&m_nameOutputFileWall);

  m_saveRateWall = 1;
  setParameter("SaveRateWall",&m_saveRateWall);

  m_alpha_function = std::string("0.");
  setParameter("Alpha",&m_alpha_function);

  m_beta_function = std::string("0.");
  setParameter("Beta",&m_beta_function);

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

  // no initialization for xref, delayed to config
  setParameter("Xref",&m_Xref);

  m_appendIter = false;
  setParameter("AppendIter",&m_appendIter);

  m_appendTime = false;
  setParameter("AppendTime",&m_appendTime);

}

//////////////////////////////////////////////////////////////////////////////

NavierStokesConsComputeAero::~NavierStokesConsComputeAero()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
NavierStokesConsComputeAero::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);
  result.push_back(&socket_states);
  result.push_back(&socket_faceNeighCell);
  result.push_back(&socket_normals);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
NavierStokesConsComputeAero::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = m_sockets.getAllSourceSockets();

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAero::setup()
{
  SafePtr<SpaceMethod> spaceMethod = getMethodData().getCollaborator<SpaceMethod>();
  SafePtr<FluctuationSplit> fs = spaceMethod.d_castTo<FluctuationSplit>();

  m_fsData = fs->getData();

  try
  {
    m_diffVar = m_fsData->getDiffusiveVar().d_castTo<NavierStokesVarSet>();
  }
  catch (FailedCastException& e)
  {
    m_diffVar = NULL;
  }

  m_updateVarSet->getModel()->resizePhysicalData(m_dataState);
  const RealVector& refValues = m_updateVarSet->getModel()->getReferencePhysicalData();
  m_cells = MeshDataStack::getActive()->getTrs("InnerCells");
  if (refValues[EulerTerm::V]   == 0. ||
      refValues[EulerTerm::RHO] == 0. ||
      refValues[EulerTerm::P]   == 0. ||
     (!is3D() && refValues[EulerTerm::VX] == 0. && refValues[EulerTerm::VY] == 0.) ||
     ( is3D() && refValues[EulerTerm::VX] == 0. && refValues[EulerTerm::VY] == 0. && refValues[EulerTerm::VZ] == 0.) )
  {
    std::string msg = std::string("Some needed reference values not set.");
    throw BadValueException (FromHere(),msg);
  }

  m_avValues.resize(m_nbeqs);

  // prepare the file for writing
  prepareOutputFileAero();
}
//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAero::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DataProcessingCom::configure(args);

  // set the dimension
  Common::SafePtr<PhysicalModel> physModel =
    PhysicalModelStack::getInstance().getEntryByNamespace(getMethodData().getNamespacePtr());

  m_dim   = physModel->getDim();
  m_nbeqs = physModel->getNbEq();

  if ( m_dim < DIM_2D || m_dim > DIM_3D ) throw BadValueException (FromHere(),"Computation of aerodynamic coefficients only valid in 2D and 3D");

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

  // set the function that computes the alpha and beta
  try
  {
    setAngleFunctions();
  }
  catch (Common::Exception& e)
  {
    CFout << e.what() << "\n";
    throw;
  }

  if (MathChecks::isZero(m_pInf))   throw BadValueException (FromHere(),"pInf, the pressure at infinity, is zero");
  if (MathChecks::isZero(m_rhoInf)) throw BadValueException (FromHere(),"rhoInf, the density at infinity, is zero");
  if (MathChecks::isZero(m_uInf))   throw BadValueException (FromHere(),"uInf, the velocity at infinity, is zero");


  CFLog(WARN, "PInf   : " << m_pInf << "\n");
  CFLog(WARN, "RhoInf : " << m_rhoInf << "\n");
  CFLog(WARN, "UInf   : " << m_uInf << "\n");

  // resize the coordinate vector
  if (m_Xref.empty())
  {
    m_Xref.resize(m_dim);
    m_Xref[XX] = 0.25;
    m_Xref[YY] = 0.0;
    if (is3D()) m_Xref[ZZ] = 0.0;
  }

  if ( m_Xref.size() != m_dim ) throw BadValueException (FromHere(),"User input Xref vector with size different than physical dimension");

  // resize vector of dimensional values
  m_dimState.resize(m_nbeqs);
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAero::setAngleFunctions()
{
  // some sanity checks
  std::vector<std::string> m_functionDef = Common::StringOps::getWords(m_alpha_function);
  if(m_functionDef.size() != 1)
   throw BadValueException (FromHere(),"NavierStokesConsComputeAero::setAngleFunctions(): alpha function wrongly defined.");

  m_alpha_func_parser.Parse(m_alpha_function, m_vars);

  if (m_alpha_func_parser.ErrorMsg() != 0)
  {
    std::string msg("ParseError in NavierStokesConsComputeAero::setAngleFunctions(): ");
    msg += std::string(m_alpha_func_parser.ErrorMsg());
    msg += " Function: " + m_alpha_function;
    msg += " Vars: "     + m_vars;
    throw Common::ParserException (FromHere(),msg);
  }

  m_functionDef = Common::StringOps::getWords(m_beta_function);
  if(m_functionDef.size() != 1)
   throw BadValueException (FromHere(),"NavierStokesConsComputeAero::setAngleFunctions(): beta function wrongly defined.");

  m_beta_func_parser.Parse(m_beta_function, m_vars);

  if (m_beta_func_parser.ErrorMsg() != 0)
  {
    std::string msg("ParseError in NavierStokesConsComputeAero::setAngleFunctions(): ");
    msg += std::string(m_beta_func_parser.ErrorMsg());
    msg += " Function: " + m_beta_function;
    msg += " Vars: "     + m_vars;
    throw Common::ParserException (FromHere(),msg);
  }
}


//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAero::executeOnTrs()
{
  CFAUTOTRACE;

  const CFuint iter = SubSystemStatusStack::getActive()->getNbIter();

  // compute wall distribution
  if(!(iter % m_saveRateWall)) { computeWall(); }

  // compute global integrated coefficients
  if(!(iter % m_saveRateAero)) { computeAero(); }
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAero::computeWall()
{
  CFAUTOTRACE;

  prepareOutputFileWall(); // file is opened here

  const CFreal R     = m_updateVarSet->getModel()->getRdim();
  const CFreal gamma = m_updateVarSet->getModel()->getGamma();
  const CFreal gammaMinus1 = gamma - 1.;

  Common::SafePtr<std::vector<CFuint> > statesIdx = getCurrentTRS()->getStatesInTrs();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  cf_assert(m_fileWall->isopen());
  ofstream& fout = m_fileWall->get();

  const CFuint rhoidx  = 0;
  const CFuint rhoUidx = 1;
  const CFuint rhoVidx = 2;
  const CFuint rhoWidx = 3;
  const CFuint rhoEidx = is3D() ? 3 : 4;

  for (CFuint iState = 0; iState < statesIdx->size(); ++iState)
  {
    const State& state = *(states[(*statesIdx)[iState]]);
    m_updateVarSet->setDimensionalValues(state,m_dimState);

    const CFreal rho  = m_dimState[rhoidx];
    const CFreal rhoU = m_dimState[rhoUidx];
    const CFreal rhoV = m_dimState[rhoVidx];
    const CFreal rhoW = is3D() ? 0. : m_dimState[rhoWidx];
    const CFreal rhoE = m_dimState[rhoEidx];

    const CFreal invRho = 1./rho;

    const CFreal rhoU2  = rhoU*rhoU;
    const CFreal rhoV2  = rhoV*rhoV;
    const CFreal rhoW2  = rhoW*rhoW;
    const CFreal U2  = rhoU2*invRho*invRho;
    const CFreal V2  = rhoV2*invRho*invRho;
    const CFreal W2  = rhoW2*invRho*invRho;

    const CFreal rhoK2  = 0.5*rho*(U2 + V2 + W2);

    // calculate pressure in dimensional values
    const CFreal P = gammaMinus1*(rhoE - rhoK2);

    // calculate pressure coefficient
    const CFreal Cp = (P-m_pInf) / (0.5*m_rhoInf*m_uInf*m_uInf);

    // calculate temperature in dimensional values
    const CFreal T = P / (rho*R);

    // calculate mach and speed of sound
    const CFreal a2 = gamma*R*T;
    const CFreal Mach = sqrt( (U2 + V2 + W2) / a2 );

    // output to File
    fout
    << (state.getCoordinates())[XX]   << " "
    << (state.getCoordinates())[YY]   << " ";

    if (is3D()) fout << (state.getCoordinates())[ZZ]   << " ";

    fout
    << Cp     << " "
    << Mach   << " "
    << P      << " "
    << T      << " "
    << rho    << "\n";
  }

  m_fileWall->close();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAero::computeAero()
{
  CFAUTOTRACE;

  // get the datahandle containing the boundary Normals
  const std::string TRSName = getCurrentTRS()->getName();
  const std::string socketName = TRSName + "-boundaryNormals";
  DataHandle<const CFreal*> boundaryNormals = m_sockets.
    getSocketSource<const CFreal*>(socketName)->getDataHandle();

  const CFreal gammaMinus1 = m_updateVarSet->getModel()->getGamma() - 1.;
  const CFreal refLength = PhysicalModelStack::getActive()->getImplementor()->getRefLength();

  CFreal tau = 0.0;
  RealVector normal(0.0, m_dim);
  RealVector Cp    (0.0, m_dim);
  RealVector Cf    (0.0, m_dim);

  // velocity variables
  RealVector u;
  RealVector v;
  RealVector w;

  // gradient velocity
  RealVector grad_u(m_dim);
  RealVector grad_v(m_dim);
  RealVector grad_w(m_dim);

  // place gradients in a vector to pass in dynamic viscosity function
  std::vector<RealVector*> vgrad;
  vgrad.push_back(&grad_u);
  vgrad.push_back(&grad_v);
  vgrad.push_back(&grad_w);

  CFreal Cmp = 0.0;
  CFreal Cmf = 0.0;
  CFreal m;
  CFreal p;

  std::vector<RealVector> vcoords(1);
  RealVector& coords = vcoords[0];
  coords.resize(m_dim);

  //  building the TRS for faces

  // builder for standard TRS GeometricEntity's that will be used for the faces
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilderFace;
  geoBuilderFace.setup();
  StdTrsGeoBuilder::GeoData& geoDataFace = geoBuilderFace.getDataGE();

  // builder for standard TRS GeometricEntity's
  Framework::GeometricEntityPool<Framework::StdTrsGeoBuilder> geoBuilderCell;
  geoBuilderCell.setup();
  StdTrsGeoBuilder::GeoData& geoDataCell = geoBuilderCell.getDataGE();

  // handle to the neighbor cell
  DataHandle<Common::Trio<CFuint, CFuint, Common::SafePtr<Framework::TopologicalRegionSet> > >
    faceNeighCell = socket_faceNeighCell.getDataHandle();

  // loop over the faces
  geoDataFace.trs = getCurrentTRS();
  const CFuint nbFaces = getCurrentTRS()->getLocalNbGeoEnts();
  cf_assert(boundaryNormals.size() == nbFaces);
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace)
  {
    /*CF_DEBUG_OBJ(iFace);*/
    geoDataFace.idx = iFace;
    GeometricEntity & currFace = *geoBuilderFace.buildGE();

    // get the face normal
    for(CFuint iDim = 0; iDim < m_dim; ++iDim)
    {
      normal[iDim] = -(boundaryNormals[iFace])[iDim];
    }
    normal.normalize();

    std::vector<State*>& states = *(currFace.getStates());
    const CFuint nb_face_states = states.size();
    cf_assert( nb_face_states >= 2);

    // Check if face is owned by one or more processors and compute the weight accordingly
    CFreal weight = 1.0;
    const CFreal weight_fraction = 1.0 / nb_face_states;
    for (CFuint s = 0; s < nb_face_states; ++s)
    {
      if (!states[s]->isParUpdatable()) { weight -= weight_fraction; }
    }

    // Compute the pressure on a face
    CFreal press = 0.0;
    m_avValues   = 0.0;
    m            = 0.0;
    coords       = 0.0;

    const CFuint rhoidx  = 0;
    const CFuint rhoUidx = 1;
    const CFuint rhoVidx = 2;
    const CFuint rhoWidx = 3;
    const CFuint rhoEidx = is3D() ? 3 : 4;

    for (CFuint s = 0; s < nb_face_states; ++s)
    {
      m_updateVarSet->setDimensionalValues(*states[s],m_dimState);
      const RealVector x = (*states[s]).getCoordinates();

      const CFreal rho  = m_dimState[rhoidx];
      const CFreal rhoU = m_dimState[rhoUidx];
      const CFreal rhoV = m_dimState[rhoVidx];
      const CFreal rhoW = is3D() ? 0. : m_dimState[rhoWidx];
      const CFreal rhoE = m_dimState[rhoEidx];

      const CFreal invRho = 1./rho;

      const CFreal rhoU2  = rhoU*rhoU;
      const CFreal rhoV2  = rhoV*rhoV;
      const CFreal rhoW2  = rhoW*rhoW;
      const CFreal U2  = rhoU2*invRho*invRho;
      const CFreal V2  = rhoV2*invRho*invRho;
      const CFreal W2  = rhoW2*invRho*invRho;

      const CFreal rhoK2  = 0.5*rho*(U2 + V2 + W2);

      p      = gammaMinus1*(rhoE - rhoK2);
      press += p;
      m     += (m_pInf - p)*((x[XX]-m_Xref[XX])*normal[YY]-(x[YY]-m_Xref[YY])*normal[XX]);

      m_avValues[rhoidx]  +=  rho;
      m_avValues[rhoUidx] +=  rhoU;
      m_avValues[rhoVidx] +=  rhoV;
      if (is3D())  m_avValues[rhoWidx] +=  rhoW;
      m_avValues[rhoEidx] +=  rhoE;

      coords += x;
    }

    press      *= weight_fraction;
    m          *= weight_fraction;
    m_avValues *= weight_fraction;
    coords     *= weight_fraction;

    const CFreal faceArea = currFace.computeVolume();

    // compute Cp, scaled by weight
    // which depends on percentage of ownership of the face within processor
    // press-m_pInf are switched to take into account the minus sign in normal
    Cp[XX] += weight * (m_pInf-press) * faceArea * normal[XX];
    Cp[YY] += weight * (m_pInf-press) * faceArea * normal[YY];
    Cp[ZZ] += weight * (m_pInf-press) * faceArea * normal[ZZ];
    Cmp    -= m*weight*faceArea ;

    // Build the neighbor cell
    const CFuint faceID = currFace.getID();
    const CFuint cellID = faceNeighCell[faceID].first;
    geoDataCell.idx = cellID;
    geoDataCell.trs = faceNeighCell[faceID].third;
    GeometricEntity& neighborCell = *geoBuilderCell.buildGE();

    // get states in the cell
    vector<Framework::State*>* states_cell_ptr = neighborCell.getStates();
    vector<Framework::State*>& states_cell = (*states_cell_ptr);
    const CFuint nbStatesInCell = states_cell.size();

    /// @todo NV: this conversion from rhou, rhov to u,v is not very good
    ///           should be done in a better way
    u.resize(nbStatesInCell);
    v.resize(nbStatesInCell);
    w.resize(nbStatesInCell);
    for (CFuint is = 0; is < nbStatesInCell ; ++is)
    {
      const State& state = (*states_cell[is]);
      const CFreal invrho = 1 / state[rhoidx];
      u[is] = state[rhoUidx] * invrho;
      v[is] = state[rhoVidx] * invrho;
      w[is] = is3D() ? 0. : state[rhoWidx] * invrho;
    }

    // Compute the shape function gradient at the average point
    const std::vector<RealMatrix> cellGradients = neighborCell.computeSolutionShapeFunctionGradients(vcoords);
    const RealMatrix& grad = cellGradients[0];

    grad_u = grad * u ;
    grad_v = grad * v ;
    grad_w = grad * w ;

    CFreal mu = 0.;
    if (m_diffVar.isNotNull()) mu = m_diffVar->getDynViscosity(m_avValues, vgrad);

    tau = mu*( normal[YY]*( grad_u[XX]*normal[XX] + grad_u[YY]*normal[YY])
             - normal[XX]*( grad_v[XX]*normal[XX] + grad_v[YY]*normal[YY]) );

    CFreal coef  = tau*weight*faceArea  ;

    Cf[XX] -= coef*normal[YY] ;
    Cf[YY] += coef*normal[XX] ;

    Cmf += (coords[YY]-m_Xref[YY])*( coef*normal[XX])
         - (coords[XX]-m_Xref[XX])*(-coef*normal[XX]);

    // release geoents
    geoBuilderFace.releaseGE();
    geoBuilderCell.releaseGE();

  } // end loop over faces

  // Compute the value of m_alpha
  m_eval[0] = SubSystemStatusStack::getActive()->getCurrentTime();
  m_alphadeg = m_alpha_func_parser.Eval(&m_eval[0]);
  m_betadeg = m_beta_func_parser.Eval(&m_eval[0]);
  
  // Transform into Radiants
  const CFreal alpharad = m_alphadeg*MathTools::MathConsts::CFrealPi()/180;
  const CFreal betarad  = m_betadeg*MathTools::MathConsts::CFrealPi()/180;

  // adimensionalize Cp
  Cp  /= (0.5 * m_rhoInf * m_uInf * m_uInf);
  Cmp /= (0.5 * m_rhoInf * m_uInf * m_uInf*refLength);
  Cf  /= (0.5 * m_rhoInf * m_uInf * m_uInf);
  Cmf /= (0.5 * m_rhoInf * m_uInf * m_uInf*refLength) ;

  const CFreal sin_a = sin(alpharad);
  const CFreal cos_a = cos(alpharad);
  const CFreal sin_b = sin(betarad);
  const CFreal cos_b = cos(betarad);

  // project Cp
  CFreal Cl   = (sin_a*(Cp[XX] + Cf[XX]) + cos_a*(Cp[YY] + Cf[YY])) / refLength;

  CFreal Cd_p = (cos_a*cos_b*Cp[XX] - sin_a*cos_b*Cp[YY] - sin_b*Cp[ZZ]) / refLength;
  CFreal Cd_f = (cos_a*cos_b*Cf[XX] - sin_a*cos_b*Cf[YY] - sin_b*Cf[ZZ]) / refLength;
  CFreal Cd   = Cd_p + Cd_f;

  CFreal Cs_p = (-cos_a*sin_b*Cp[XX] + sin_a*sin_b*Cp[YY] - cos_b*Cp[ZZ]) / refLength;
  CFreal Cs_f = (-cos_a*sin_b*Cf[XX] + sin_a*sin_b*Cf[YY] - cos_b*Cf[ZZ]) / refLength;
  CFreal Cs   = Cs_p + Cs_f;

  CFreal actual_momentum_p = Cmp / refLength;
  CFreal actual_momentum_f = Cmf / refLength;
  CFreal actual_momentum = Cmp + Cmf;

  /// @TODO NV: has been copied from RMSJouleHeatSource.cxx
  ///           where this solution was said to be temporary
  CFreal total_Cl   = 0.0;
  CFreal total_Cd_p = 0.0;
  CFreal total_Cd_f = 0.0;
  CFreal total_Cd   = 0.0;
  CFreal total_Cs   = 0.0;
  CFreal total_momentum_p = 0.0;
  CFreal total_momentum_f = 0.0;
  CFreal total_momentum = 0.0;
  
  const std::string nsp = getMethodData().getNamespace();

  if (PE::GetPE().GetProcessorCount(nsp) > 1) {

#ifdef CF_HAVE_MPI
  MPI_Allreduce(&Cl,   &total_Cl,   1, MPI_DOUBLE, MPI_SUM, PE::GetPE().GetCommunicator(nsp));
  MPI_Allreduce(&Cd_p, &total_Cd_p, 1, MPI_DOUBLE, MPI_SUM, PE::GetPE().GetCommunicator(nsp));
  MPI_Allreduce(&Cd_f, &total_Cd_f, 1, MPI_DOUBLE, MPI_SUM, PE::GetPE().GetCommunicator(nsp));
  MPI_Allreduce(&Cd,   &total_Cd,   1, MPI_DOUBLE, MPI_SUM, PE::GetPE().GetCommunicator(nsp));
  MPI_Allreduce(&Cs,   &total_Cs,   1, MPI_DOUBLE, MPI_SUM, PE::GetPE().GetCommunicator(nsp));
  MPI_Allreduce(&actual_momentum_p, &total_momentum_p, 1, MPI_DOUBLE, MPI_SUM, PE::GetPE().GetCommunicator(nsp));
  MPI_Allreduce(&actual_momentum_f, &total_momentum_f, 1, MPI_DOUBLE, MPI_SUM, PE::GetPE().GetCommunicator(nsp));
  MPI_Allreduce(&actual_momentum,   &total_momentum,   1, MPI_DOUBLE, MPI_SUM, PE::GetPE().GetCommunicator(nsp));
#endif
  
  }
  else
  {

    total_Cl   = Cl;
    total_Cd_p = Cd_p;
    total_Cd_f = Cd_f;
    total_Cd   = Cd;
    total_Cs   = Cs;
    total_momentum_p = actual_momentum_p;
    total_momentum_f = actual_momentum_f;
    total_momentum   = actual_momentum;

  }

  // Output to file the coefficients
  updateOutputFileAero();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAero::unsetup()
{
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAero::prepareOutputFileWall()
{


  cf_assert (!m_fileWall->isopen());

  boost::filesystem::path file = Environment::DirPaths::getInstance().getResultsDir() /
boost::filesystem::path ( m_nameOutputFileWall );
  file = PathAppender::getInstance().appendAllInfo(file, m_appendIter, m_appendTime );

  ofstream& fout = m_fileWall->open(file);

  fout << "TITLE  =  Values at the Wall" << "\n";
  fout << "VARIABLES = x y Cp Mach Pressure Temperature Density" << "\n";
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAero::prepareOutputFileAero()
{
  boost::filesystem::path fpath = m_nameOutputFileAero;
  fpath = PathAppender::getInstance().appendAllInfo(fpath, false, false);

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& convergenceFile = fhandle->open(fpath);

  convergenceFile <<"TITLE  =  Aerodynamic Coefficients"  << "\n";
  convergenceFile << "VARIABLES = Iter PhysTime Alpha CL CDp CDf CD Cs CMp CMf CM" << "\n";
  convergenceFile.close();

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void NavierStokesConsComputeAero::updateOutputFileAero()
{
  boost::filesystem::path fpath (m_nameOutputFileAero);
  fpath = PathAppender::getInstance().appendAllInfo( fpath, false, false);

  SelfRegistPtr<Environment::FileHandlerOutput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerOutput>::getInstance().create();
  ofstream& convergenceFile = fhandle->open(fpath, ios_base::app);

  Common::SafePtr<SubSystemStatus> subSysStatus =
    SubSystemStatusStack::getActive();

  convergenceFile
  << subSysStatus->getNbIter()       << " "
  << subSysStatus->getCurrentTime()  << " "
  << m_alphadeg                      << " "
  << total_Cl                      << " "
  << total_Cd_p                    << " "
  << total_Cd_f                    << " "
  << total_Cd                      << " "
  << total_Cs                      << " "
  << total_momentum_p                << " "
  << total_momentum_f                << " "
  << total_momentum                  << "\n";

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NavierStokes

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

