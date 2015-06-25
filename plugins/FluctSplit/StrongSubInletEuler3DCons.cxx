#include <numeric>
#include <algorithm>

#include "StrongSubInletEuler3DCons.hh"
#include "FluctSplit/CreateBoundaryNodalNormals.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "FluctSplit/FluctSplitSpaceTimeNavierStokes.hh"
#include "Framework/MeshData.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "MathTools/MatrixInverter.hh"
#include "Environment/ObjectProvider.hh"
#include "FluctSplit/SpaceTime_Splitter.hh"
#include "Framework/CFL.hh"
#include "MathTools/RealMatrix.hh"
#include "Framework/SubSystemStatus.hh"
#include "Common/BadValueException.hh"
#include "Common/ParserException.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StrongSubInletEuler3DCons, FluctuationSplitData, FluctSplitSpaceTimeNavierStokesModule> 
strongSubFunctionInletLinEuler3DConsProvider("StrongSubInletEuler3DCons");

//////////////////////////////////////////////////////////////////////////////

StrongSubInletEuler3DCons::StrongSubInletEuler3DCons(const std::string& name) :
  FluctuationSplitCom(name),
  socket_normals("normals"),
  socket_faceNeighCell("faceNeighCell"),
  socket_rhs("rhs"),
  socket_states("states"),
  socket_isUpdated("isUpdated"),
  socket_updateCoeff("updateCoeff"),
  socket_nodes("nodes"),
  socket_isBState("isBState"), 
  _varSet(),
  m_nbEqs(),
  _r1(),
  _r2(),
  _r3(),
  _r4()
{
  addConfigOptionsTo(this);

  m_vars_inflow.resize(0);

  m_vars_inflow = std::vector<std::string>();
  setParameter("Vars",&m_vars_inflow);

  m_function_inflow = vector<std::string>();
  setParameter("InFlow",&m_function_inflow);
  
  m_useBlasius = false;
  setParameter("UseBlasiusInflow",&m_useBlasius);
  
  m_ReferenceLength = 1.0;
  setParameter("ReferenceLength",&m_ReferenceLength);
  
  m_UpstreamPosition = 1.0;
  setParameter("UpstreamPosition",&m_UpstreamPosition);
  
  m_ReynoldsNumber = 1.0;
  setParameter("ReynoldsNumber",&m_ReynoldsNumber);
  
  m_MachNumber = 1.0;
  setParameter("MachNumber",&m_MachNumber);
  
  m_gamma = 1.4;
  setParameter("Gamma",&m_gamma);
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubInletEuler3DCons::defineConfigOptions(Config::OptionList&
options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Variable names.");
  options.addConfigOption< std::vector<std::string> >("InFlow","Function defining the incoming flow.");
  options.addConfigOption< bool >("UseBlasiusInflow","Generates a Blasius inflow profile based on the distance to the wall, Reynolds number and upstream position");
  options.addConfigOption< CFreal >("ReferenceLength","Reference Length of the system [m].");
  options.addConfigOption< CFreal >("UpstreamPosition","Position [m] aft of the leading edge at which the Blasius profile should be inscribed.");
  options.addConfigOption< CFreal >("ReynoldsNumber","Reynolds number based on reference length, freestream velocity and viscosity.");
  options.addConfigOption< CFreal >("MachNumber","Mach number based on reference length, freestream velocity and viscosity.");
  options.addConfigOption< CFreal >("Gamma","Adiabatic exponent");
}

//////////////////////////////////////////////////////////////////////////////

StrongSubInletEuler3DCons::~StrongSubInletEuler3DCons()
{
  if (isSetup()) unsetup();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StrongSubInletEuler3DCons::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_normals);
  result.push_back(&socket_faceNeighCell);
  result.push_back(&socket_rhs);
  result.push_back(&socket_states);
  result.push_back(&socket_isUpdated);
  result.push_back(&socket_updateCoeff);
  result.push_back(&socket_nodes);
  result.push_back(&socket_isBState);  
  

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubInletEuler3DCons::setup()
{
  FluctuationSplitCom::setup();

  m_nbEqs = PhysicalModelStack::getActive()->getNbEq();

// create boundary nodal normals, which pointing outwards
  _bcNormals.resize(getTrsList().size());
  _varSet->setup();
  
  _r1.resize(PhysicalModelStack::getActive()->getNbEq());
  _r2.resize(PhysicalModelStack::getActive()->getNbEq());
  _r3.resize(PhysicalModelStack::getActive()->getNbEq());
  _r4.resize(PhysicalModelStack::getActive()->getNbEq());

  CreateBoundaryNodalNormals obj(getMethodData().getStdTrsGeoBuilder());
  obj.setDataSockets(socket_normals,socket_faceNeighCell);
  obj.create(getTrsList(), _bcNormals);


// handling the inflow disturbances
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  inflow.resize(states.size());

  const CFuint nb_funs = m_function_inflow.size();

  for (CFuint i = 0; i < states.size(); ++i)
  {
    inflow[i].resize(m_nbEqs);
    for (CFuint j = 0; j < m_nbEqs; ++j)
      inflow[i][j] = 0.0;
  }

/// create the cell connectivity

  // find the inner trs
  vector< SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  SafePtr<TopologicalRegionSet> innerTrs=0;
  for (CFuint i = 0; i < trs.size(); ++i) {
    if (trs[i]->hasTag("inner")) {
      innerTrs = trs[i];
      break;
    }
  }
  if (innerTrs==0) Common::NoSuchValueException(FromHere(),"Trs with tag 'inner' not found.");
  // set up numbers
  DataHandle< Node*,GLOBAL > innerNodes = socket_nodes.getDataHandle();
  const CFuint innerNbNodes = innerNodes.size();
  CFVec<CFint> bndNodeGlobal2Local(innerNbNodes);
  const CFuint innerNbGeos=innerTrs->getLocalNbGeoEnts();
  const Common::SafePtr< Common::ConnectivityTable<CFuint> > innerGeo2Node=innerTrs->getGeo2NodesConn();


  // set the boundary trs
  SafePtr<TopologicalRegionSet> bndTrs = getCurrentTRS();
   // getting nodes and preliminary settings
  Common::SafePtr< std::vector< CFuint > > bndNodes = bndTrs->getNodesInTrs();
  const CFuint bndNbNodes= bndTrs->getNbNodesInTrs();
  bndNod2Elm.resize(bndNbNodes);
  bndNodeGlobal2Local=-1;
  // cycle all the nodes in the TRS to set isBndNode flag
  CFuint bndCheckNbNode=0;
  for (std::vector< CFuint >::iterator itd = bndNodes->begin(); itd != bndNodes->end(); ++itd)
     bndNodeGlobal2Local[*itd]=bndCheckNbNode++;

  // build boundary node -> innercells elem connectivity
  for (CFuint i=0; i<innerNbGeos; i++) {
    const CFuint innerGeoLocalID=innerTrs->getLocalGeoID(i);
    const CFuint innerNbGeoNodes=innerGeo2Node->nbCols(i);
    for (CFuint j=0; j<innerNbGeoNodes; j++)
      if (bndNodeGlobal2Local[(*innerGeo2Node)(i,j)]!=-1)
        bndNod2Elm[bndNodeGlobal2Local[(*innerGeo2Node)(i,j)]].push_back(innerGeoLocalID);
  }
  
  
  
  Common::SafePtr< vector<CFuint> > const statesIdx = getCurrentTRS()->getStatesInTrs();
  blasius_inflow.resize(states.size());
  blasius_energy_inflow.resize(states.size());
  
// Evaluate the flow variables from the input formula/values
  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
   const CFuint stateID = (*statesIdx)[iState]; 
   
    if(m_useBlasius)
    {
      /// Calculate the Blasius velocity profile for the different positions
      /// Uses RK4 
	
      Node& coord = states[stateID]->getCoordinates();
      if(coord[1]<0.0) std::cout<<"ERROR, y<0"<<std::endl;
      CFreal y1		=	0.0;
      CFreal y2		=	0.0;
      CFreal y3		=	0.4696;
      
      CFreal eta 		=       0.0;
      CFreal h	=	0.01;	
      CFreal eta_end	=	coord[1]/sqrt(2.0/m_ReynoldsNumber*m_UpstreamPosition/m_ReferenceLength);
      CFreal y1temp,y2temp,y3temp;
      CFreal k11,k12,k13;
      CFreal k21,k22,k23;
      CFreal k31,k32,k33;
      CFreal k41,k42,k43;
      
      if(eta_end>10.0) y2=1.00;
      else
      {
	while(eta<eta_end)
	{
	  k11=h*y2; k12=h*y3; k13	=h*(-y3*y1);
	  y1temp = y1+0.5*k11; y2temp = y2+0.5*k12; y3temp= y3+0.5*k13;
	  
	  k21=h*y2temp; k22=h*y3temp; k23 = h*(-y1temp*y3temp);
	  y1temp= y1+0.5*k21; y2temp = y2+0.5*k22; y3temp = y3+0.5*k23;
	  
	  k31=h*y2temp; k32=h*y3temp; k33=h*(-y1temp*y3temp);
	  y1temp= y1+k31; y2temp = y2+k32; y3temp = y3+k33;
	  
	  k41=h*y2temp; k42=h*y3temp; k43=h*(-y1temp*y3temp);
	  
	  y1=y1+(k11+2.0*k21+2.0*k31+k41)/6.0;
	  y2=y2+(k12+2.0*k22+2.0*k32+k42)/6.0;
	  y3=y3+(k13+2.0*k23+2.0*k33+k43)/6.0;
	  eta = eta + h;
	}
      }
      blasius_inflow[stateID]=y2;
      blasius_energy_inflow[stateID]=1.0/m_gamma/(m_gamma-1.0)/m_MachNumber/m_MachNumber;
      std::cout<<"y="<<coord[1]<<" eta="<<eta-h<<" f'(eta)="<<y2<<" E-k="<<blasius_energy_inflow[stateID]<<std::endl;
   }
}
/// /*-----------------------------------------------------------------------*/
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubInletEuler3DCons::executeOnTrs()
{
  vector<RealVector>* bcNormalsInTrs = &_bcNormals[getCurrentTrsID()];
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<bool> isUpdated = socket_isUpdated.getDataHandle();
  DataHandle<CFreal> updateCoeff = socket_updateCoeff.getDataHandle();
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();
  
  DataHandle< InwardNormalsData*> m_normals = socket_normals.getDataHandle();

  Common::SafePtr< vector<CFuint> > const statesIdx = getCurrentTRS()->getStatesInTrs();
  
  const RealVector& linearData = _varSet->getModel()->getPhysicalData();
   
  const CFuint dim  = PhysicalModelStack::getActive()->getDim();
  const CFuint TT   = dim;
  const CFreal time =  SubSystemStatusStack::getActive()->getCurrentTimeDim();
  m_var_values[TT] = time;
  
  const CFreal gamma = _varSet->getModel()->getGamma();


// Evaluate the flow variables from the input formula/values
  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
   const CFuint stateID = (*statesIdx)[iState]; 
   Node& coord = states[stateID]->getCoordinates();
   for (CFuint iCoor = 0; iCoor < dim; ++iCoor) {
    m_var_values[iCoor] = coord[iCoor];
  }
    m_function_parser_inflow.evaluate(m_var_values, inflow[stateID]);
    
    /// Apply Blasius boundary layer profile that has been calculated in setup
    if(m_useBlasius)
    {
      //Node& coord = states[stateID]->getCoordinates();
      inflow[stateID][1]=inflow[stateID][1]*blasius_inflow[stateID];
      inflow[stateID][4]=inflow[stateID][0]*blasius_energy_inflow[stateID]+0.5/inflow[stateID][0]*inflow[stateID][1]*inflow[stateID][1];
    }

  }
      
  for (CFuint iState = 0; iState < statesIdx->size(); ++iState) {
    const CFuint stateID = (*statesIdx)[iState];
    State *state = states[stateID];
    const RealVector* bcNormal = &(*bcNormalsInTrs)[stateID];
    
//     get the boundary elements containing the state
    vector <CFuint> cells = bndNod2Elm[iState];
    CFuint nbElemsOfState = cells.size();

//     RealVector bcNormalState =  bndNodNorm[iState];

    // oriented to the other direction
//     CFreal nx = -(bcNormalState)[0];
//     CFreal ny = -(bcNormalState)[1];
//     CFreal nz = -(bcNormalState)[2];
    
    const CFreal nx = (*bcNormal)[0];
    const CFreal ny = (*bcNormal)[1];
    const CFreal nz = (*bcNormal)[2];    
    
//     cout << nx << " " << ny << " " << nz << "\n";
     
    (*state)[0] = inflow[stateID][0];
    (*state)[1] = inflow[stateID][1];
    (*state)[2] = inflow[stateID][2];
    (*state)[3] = inflow[stateID][3];
    (*state)[4] = inflow[stateID][4];   

    if (!isUpdated[stateID]) {
      const CFreal updateCoeffValue = updateCoeff[stateID];
      if (std::abs(updateCoeffValue) > MathTools::MathConsts::CFrealEps()) {
	
  _varSet->setEigenVect1(_r1, *state, *bcNormal);
  _varSet->setEigenVect2(_r2, *state, *bcNormal);
  _varSet->setEigenVect3(_r3, *state, *bcNormal);
  _varSet->setEigenVect4(_r4, *state, *bcNormal);
	
  CFreal *const rhsStart = &rhs(stateID, 0, m_nbEqs);

  CFreal drho = rhsStart[0];
  CFreal drho0u = rhsStart[1];
  CFreal drho0v = rhsStart[2];
  CFreal drho0w = rhsStart[3];
  CFreal drhoE = rhsStart[4];


  const CFreal avRho = linearData[EulerTerm::RHO];
  const CFreal avU   = linearData[EulerTerm::VX];
  const CFreal avV   = linearData[EulerTerm::VY];
  const CFreal avW   = linearData[EulerTerm::VZ];  
  const CFreal avH   = linearData[EulerTerm::H];
  const CFreal avA   = linearData[EulerTerm::A];
  const CFreal ovAvA = 1./avA;
  const CFreal ovAvA2 = ovAvA*ovAvA;
  const CFreal gammaMinus1 = gamma - 1.;
  const CFreal um = avU*nx + avV*ny + avW*nz;
  const CFreal ra = 0.5*avRho*ovAvA;
  const CFreal k = 0.5*(avU*avU + avV*avV +  avW*avW);
  const CFreal k1 = 1.-k*gammaMinus1*ovAvA2;
  const CFreal k2 = -gammaMinus1*ovAvA2;
  const CFreal ovAvRho = 1./avRho;
  const CFreal uDivA = gammaMinus1*avU*ovAvA;
  const CFreal vDivA = gammaMinus1*avV*ovAvA;
  const CFreal wDivA = gammaMinus1*avW*ovAvA;  
  const CFreal ovAvRhoA = ovAvRho*ovAvA;
	
  const CFreal beta1 = (nx*k1-ovAvRho*(nz*avV-ny*avW))*drho + nx*uDivA*ovAvA*drho0u +  (nx*vDivA*ovAvA+nz*ovAvRho)*drho0v + (nx*wDivA*ovAvA-ny*ovAvRho)*drho0w + k2*nx*drhoE;
  const CFreal beta2 = (ny*k1-ovAvRho*(nx*avW-nz*avU))*drho + (ny*uDivA*ovAvA-nz*ovAvRho)*drho0u +  (ny*vDivA*ovAvA)*drho0v + (ny*wDivA*ovAvA+nx*ovAvRho)*drho0w + k2*ny*drhoE;
  const CFreal beta3 = (nz*k1-ovAvRho*(ny*avU-nx*avV))*drho + (nz*uDivA*ovAvA+ny*ovAvRho)*drho0u +  (nz*vDivA*ovAvA-nx*ovAvRho)*drho0v + (nz*wDivA*ovAvA)*drho0w + k2*nz*drhoE;
  const CFreal beta4 = avA*ovAvRho*(-k*gammaMinus1*ovAvA2-um*ovAvA)*drho + ovAvRho*(nx-uDivA)*drho0u + ovAvRho*(ny-vDivA)*drho0v + ovAvRho*(nz-wDivA)*drho0w + gammaMinus1* ovAvRhoA*drhoE;
	
	
        _r1 *= -beta1;
        _r2 *= -beta2;
        _r3 *= -beta3;
        _r4 *= -beta4;
        _r4 += _r3 += _r2 += _r1; // This is actually the whole boundary correction term

        for (CFuint iEq = 0; iEq < m_nbEqs; ++iEq) {
         rhs(stateID, iEq, m_nbEqs) += _r4[iEq];
        }
      }
      isUpdated[stateID] = true; // flagging is important!!!!!
    }
  }

  
}

//////////////////////////////////////////////////////////////////////////////

void StrongSubInletEuler3DCons::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);

  std::string name = getMethodData().getNamespace();
  Common::SafePtr<Namespace> nsp = NamespaceSwitcher::getInstance
    (SubSystemStatusStack::getCurrentName()).getNamespace(name);
  Common::SafePtr<PhysicalModel> physModel = PhysicalModelStack::getInstance().getEntryByNamespace(nsp);

  std::string varSetName = "Euler3DCons";
  _varSet.reset((Environment::Factory<ConvectiveVarSet>::getInstance().getProvider(varSetName)->
    create(physModel->getImplementor()->getConvectiveTerm())).d_castTo<Physics::NavierStokes::Euler3DCons>());

  cf_assert(_varSet.isNotNull());

  m_var_values.resize(m_vars_inflow.size());

  if(m_function_inflow.empty())
     throw BadValueException(FromHere(),"StrongSubInletEuler3DCons::setFuntion(): no incoming flow function provided.");

  // configure the expression for the mean flow
  m_function_parser_inflow.setFunctions(m_function_inflow);
  m_function_parser_inflow.setVariables(m_vars_inflow);

  try
  {
    m_function_parser_inflow.parse();
  }
  catch (Common::ParserException& e)
  {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }

  

}

//////////////////////////////////////////////////////////////////////////////
void StrongSubInletEuler3DCons::unsetup()
{
   inflow.resize(0);
   _bcNormals.resize(0);
  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
