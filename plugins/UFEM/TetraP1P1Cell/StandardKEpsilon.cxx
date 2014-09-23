#include "Environment/ObjectProvider.hh"
#include "Common/BadValueException.hh"

#include "UFEM/TetraP1P1Cell/StandardKEpsilon.hh"
#include "UFEM/TetraP1P1Cell/UFEMTetraP1P1Cell.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace UFEM {
    namespace TetraP1P1Cell {

using namespace std;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider < StandardKEpsilon,
                              UFEMTerm,
                              UFEMTetraP1P1CellPlugin,
                              UFEMTerm::NARGS >
aTetraP1P1Cell_StandardKEpsilon_Provider ( "TetraP1P1Cell_StandardKEpsilon" );

//////////////////////////////////////////////////////////////////////////////

void StandardKEpsilon::defineConfigOptions ( Config::OptionList& options )
{
  CFAUTOTRACE;

  options.addConfigOption< CFreal >("Cmu",        "Cmu");
  options.addConfigOption< CFreal >("Ceps1",      "Cepsilon1");
  options.addConfigOption< CFreal >("Ceps2",      "Cepsilon2");
  options.addConfigOption< CFreal >("SigmaK",     "Turbulent Prandtl Number for K");
  options.addConfigOption< CFreal >("SigmaEPS",   "Turbulent Prandtl Number for EPSILON");
  options.addConfigOption< CFreal >("MuLam",      "Laminar viscosity");
  options.addConfigOption< CFreal >("TopLimit",   "Turb/Laminar viscosity ratio, above the turbulent viscosity is cropped.");
  options.addConfigOption< CFreal >("BottomLimit","Turb/Laminar viscosity ratio, below the turbulent viscosity is cropped.");
  options.addConfigOption< CFreal >("TopPrec",    "Max value limiter of k and epsilon.");
  options.addConfigOption< CFreal >("BottomPrec", "Min value limiter of k and epsilon.");
}


//////////////////////////////////////////////////////////////////////////////

StandardKEpsilon::StandardKEpsilon ( const std::string& name, Common::SafePtr<ElemProps> props  ) :
  UFEMTerm ( name, props ),
  m_cellprops ( *props.d_castTo<TetraP1P1Cell::CellProps>() ),
  socket_interStates("interStates"),
  socket_wallNearestDistance("wallNearestDistance"),
  socket_wallNearestSegment("wallNearestSegment"),
  socket_wallNearestVelocityGradient("wallNearestVelocityGradient"),
  socket_connBndFace2InnerCell("connBndFace2InnerCell")
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  //Constants Default Values
  m_Cmu         = 0.09;
  m_Ceps1       = 1.44;
  m_Ceps2       = 1.92;
  m_SigmaK      = 1.;
  m_SigmaEPS    = 1.3;
  m_MuLam       = 1.;
  m_TopLimit    = 1.e+5;
  m_BottomLimit = 1.e-5;
  m_TopPrec     = 1.e+30;
  m_BottomPrec  = 1.e-30;

  setParameter( "Cmu",          &m_Cmu );
  setParameter( "Ceps1",        &m_Ceps1 );
  setParameter( "Ceps2",        &m_Ceps2 );
  setParameter( "SigmaK",       &m_SigmaK );
  setParameter( "SigmaEPS",     &m_SigmaEPS );
  setParameter( "MuLam",        &m_MuLam );
  setParameter( "TopLimit",     &m_TopLimit );
  setParameter( "BottomLimit",  &m_BottomLimit );
  setParameter( "TopPrec",      &m_TopPrec );
  setParameter( "BottomPrec",   &m_BottomPrec );
}

//////////////////////////////////////////////////////////////////////////////

StandardKEpsilon::~StandardKEpsilon() 
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void StandardKEpsilon::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  UFEMTerm::configure ( args );

  CFLog(INFO, getClassName() << ": Cmu: "                                      << m_Cmu          << "\n");
  CFLog(INFO, getClassName() << ": Cepsilon1: "                                << m_Ceps1        << "\n");
  CFLog(INFO, getClassName() << ": Cepsilon2: "                                << m_Ceps2        << "\n");
  CFLog(INFO, getClassName() << ": Turbulent Prandtl Number for K: "           << m_SigmaK       << "\n");
  CFLog(INFO, getClassName() << ": Turbulent Prandtl Number for EPSILON: "     << m_SigmaEPS     << "\n");
  CFLog(INFO, getClassName() << ": Laminar viscosity: "                        << m_MuLam        << "\n");
  CFLog(INFO, getClassName() << ": TopLimit (viscosity limit by ratio): "      << m_TopLimit     << "\n");
  CFLog(INFO, getClassName() << ": BottomLimit (viscosity limit by ratio): "   << m_BottomLimit  << "\n");
  CFLog(INFO, getClassName() << ": TopPrec (max of k and epsilon): "           << m_TopPrec      << "\n");
  CFLog(INFO, getClassName() << ": BottomPrec (min of k and epsilon): "        << m_BottomPrec   << "\n");

}

//////////////////////////////////////////////////////////////////////////////

void StandardKEpsilon::setup ()
{
  CFAUTOTRACE;
  UFEMTerm::setup ();

  socket_interStates.setParentNamespace ( getMethodData().getNamespace() );
  socket_wallNearestDistance.setParentNamespace ( getMethodData().getNamespace() );
  socket_wallNearestSegment.setParentNamespace ( getMethodData().getNamespace() );
  socket_wallNearestVelocityGradient.setParentNamespace ( getMethodData().getNamespace() );
  socket_connBndFace2InnerCell.setParentNamespace ( getMethodData().getNamespace() );

}

//////////////////////////////////////////////////////////////////////////////

void StandardKEpsilon::unsetup ()
{
  CFAUTOTRACE;
  UFEMTerm::unsetup ();
}

//////////////////////////////////////////////////////////////////////////////

void StandardKEpsilon::pre ()
{
  CFAUTOTRACE;
  m_KBottomLimitCounter=0;
  m_KTopLimitCounter=0;
  m_EBottomLimitCounter=0;
  m_ETopLimitCounter=0;
  
  
  // calc of wall shear stress

  // 1. necessary stuff
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  DataHandle< CFuint > wallNearestSegment = socket_wallNearestSegment.getDataHandle();
  DataHandle< CFreal > wallNearestDistance = socket_wallNearestDistance.getDataHandle();
  DataHandle< CFreal > wallNearestVelocityGradient = socket_wallNearestVelocityGradient.getDataHandle();
  DataHandle< std::vector<CFuint> > connBndFace2InnerCell = socket_connBndFace2InnerCell.getDataHandle();
  DataHandle<Framework::State*> interStates = socket_interStates.getDataHandle();
  SafePtr< GeometricEntityPool< StdTrsGeoBuilder > > geoBuilder = getMethodData().getStdTrsGeoBuilder();
  StdTrsGeoBuilder::GeoData& geoData = geoBuilder->getDataGE();

  // 2. finding trs named as "wall" and the trs of innercells
  std::vector<Common::SafePtr<TopologicalRegionSet> > trs = MeshDataStack::getActive()->getTrsList();
  SafePtr< TopologicalRegionSet > innerTrs=0;
  CFuint wallNearestTrs=-1;
  for (CFuint i = 0; i < trs.size(); ++i) {
    if (trs[i]->hasTag("inner")) {
      innerTrs = trs[i];
    }
    if (trs[i]->getName() == "wall" ) {
      wallNearestTrs=i;
    }
  }
  if (innerTrs==0) Common::BadValueException(FromHere(),"No TRS found with tag 'inner'.");
  geoData.trs = innerTrs;

  // 3. looping on bnd elements
  CFuint nbBndGeoEnt=wallNearestVelocityGradient.size();
  for (CFuint i = 0; i < nbBndGeoEnt; ++i) {
    CFuint iCell= connBndFace2InnerCell[wallNearestTrs][i];
    geoData.idx = iCell;
    const GeometricEntity& cell = *geoBuilder->buildGE();
    const vector< State* >& states = cell.getStates();
    RealVector estate(nbEqs);
    estate = 0.;
    for (CFuint iState=0; iState<4; ++iState) {
      State& state = *states[iState];
      State& interState = *interStates[state.getLocalID()];
      for (CFuint iEq=0; iEq<nbEqs; ++iEq) {
        estate[iEq] += interState[iEq];
      }
    }
    for (CFuint iEq=0; iEq<nbEqs; ++iEq) estate[iEq] *= 1./4.;
    CFreal Umag=sqrt(estate[1]*estate[1]+estate[2]*estate[2]+estate[3]*estate[3]);
    wallNearestVelocityGradient[i]=Umag/wallNearestDistance[iCell];

/*
const std::vector<Node*>& nodes = cell.getNodes();
Node& node0 = *nodes[0];
Node& node1 = *nodes[1];
Node& node2 = *nodes[2];
Node& node3 = *nodes[3];
CFLog(INFO,"y coord=" << (node0[1]+node1[1]+node2[1]+node3[1])/4. << " dist=" << wallNearestDistance[iCell] << " wallgrad=" << wallNearestVelocityGradient[i] << "\n");
*/

    geoBuilder->releaseGE();
  }

  
}

//////////////////////////////////////////////////////////////////////////////

void StandardKEpsilon::compute ( const Framework::GeometricEntity& cell , AssemblyData& adata )
{
  CFAUTOTRACE;

#define kronecker(i,j) ((i)==(j)?1.:0.)

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFreal RhoElm = getMethodData().getRhoElm();

  RealVector estate(nbEqs);
  estate = 0.;

  const vector< State* >& states = cell.getStates();
  DataHandle<Framework::State*> interStates = socket_interStates.getDataHandle();
  DataHandle< CFuint > wallNearestSegment = socket_wallNearestSegment.getDataHandle();
  DataHandle< CFreal > wallNearestDistance = socket_wallNearestDistance.getDataHandle();
  DataHandle< CFreal > wallNearestVelocityGradient = socket_wallNearestVelocityGradient.getDataHandle();

  CFreal dudx=0.,dudy=0.,dudz=0.,dvdx=0.,dvdy=0.,dvdz=0.,dwdx=0.,dwdy=0.,dwdz=0.,dKdx=0.,dKdy=0.,dKdz=0.,dEPSdx=0.,dEPSdy=0.,dEPSdz=0.;

  const CFreal toplimit = m_MuLam*m_TopLimit/getMethodData().getRhoElm();
  const CFreal bottomlimit = m_MuLam*m_BottomLimit/getMethodData().getRhoElm();

  double invvol=1./(3.*m_cellprops.getCellData().vol);
  const CFreal* nx=m_cellprops.getCellData().nx;
  const CFreal* ny=m_cellprops.getCellData().ny;
  const CFreal* nz=m_cellprops.getCellData().nz;
  for (CFuint iState=0; iState<4; ++iState) {
    State& state = *states[iState];
    State& interState = *interStates[state.getLocalID()];

    const CFreal mut=RhoElm*m_Cmu*interState[4]*interState[4]/interState[5];

    if ((mut<bottomlimit)||(interState[4]<m_BottomPrec)||(interState[5]<m_BottomPrec)) {
      if (interState[4]<m_BottomPrec) interState[4]=m_BottomPrec;
      if (interState[5]<m_BottomPrec) interState[5]=m_BottomPrec;
      const CFreal klime=sqrt(bottomlimit*interState[5]/m_Cmu);
      if (interState[4]<klime) {
        // no change interState[4]
        interState[5]=m_Cmu*interState[4]*interState[4]/bottomlimit;
        //cout << "k ";
        m_KBottomLimitCounter++;
      } else {
        interState[4]=klime;
        // no change interState[5]
        //cout << "e ";
        m_EBottomLimitCounter++;
      }
    }

    if ((mut>toplimit)||(interState[4]>m_TopPrec)||(interState[5]>m_TopPrec)){
      if (interState[4]>m_TopPrec) interState[4]=m_TopPrec;
      if (interState[5]>m_TopPrec) interState[5]=m_TopPrec;
      const CFreal klime=sqrt(toplimit*interState[5]/m_Cmu);
      if (interState[4]<klime) {
        // no change interState[4]
        interState[5]=m_Cmu*interState[4]*interState[4]/bottomlimit;
        //cout << "K ";
        m_KTopLimitCounter++;
      } else {
        interState[4]=klime;
        // no change interState[5]
        //cout << "E ";
        m_ETopLimitCounter++;
      }
    }

    dudx      +=interState[1]*nx[iState];
    dudy      +=interState[1]*ny[iState];
    dudz      +=interState[1]*nz[iState];
    dvdx      +=interState[2]*nx[iState];
    dvdy      +=interState[2]*ny[iState];
    dvdz      +=interState[2]*nz[iState];
    dwdx      +=interState[3]*nx[iState];
    dwdy      +=interState[3]*ny[iState];
    dwdz      +=interState[3]*nz[iState];
    dKdx      +=interState[4]*nx[iState];
    dKdy      +=interState[4]*ny[iState];
    dKdz      +=interState[4]*nz[iState];
    dEPSdx    +=interState[5]*nx[iState];
    dEPSdy    +=interState[5]*ny[iState];
    dEPSdz    +=interState[5]*nz[iState];
    for (CFuint iEq=0; iEq<nbEqs; ++iEq) {
      estate[iEq] += interState[iEq];
    }
  }
  dudx   *= invvol;
  dudy   *= invvol;
  dudz   *= invvol;
  dvdx   *= invvol;
  dvdy   *= invvol;
  dvdz   *= invvol;
  dwdx   *= invvol;
  dwdy   *= invvol;
  dwdz   *= invvol;
  dKdx   *= invvol;
  dKdy   *= invvol;
  dKdz   *= invvol;
  dEPSdx *= invvol;
  dEPSdy *= invvol;
  dEPSdz *= invvol;
  for (CFuint iEq=0; iEq<nbEqs; ++iEq) estate[iEq] *= 1./4.;

  // assign understandible names for it
  const CFreal K           = estate[3];
  const CFreal epsilon     = estate[4];

  CFuint cid=cell.getID();

  // low-Re wall damping
  
  // Chien
  m_Ceps1  = 1.35;
  m_Ceps2  = 1.8;
  CFreal Rt = K*K/(m_MuLam*epsilon);  // !!! some articles say to multiply with Re_infinity
  CFreal tauwall = m_MuLam*wallNearestVelocityGradient[wallNearestSegment[cid]];
  CFreal ustar = sqrt(tauwall/RhoElm);
  CFreal yplus = RhoElm*ustar*wallNearestDistance[cid]/m_MuLam;
  CFreal fmu = 1-exp(-0.0115*yplus);
  CFreal feps1 = 1.;
  CFreal feps2 = 1-0.22*exp(-(Rt/6.)*(Rt/6.));

/*
  // Lam & Bremhorst  
  CFreal Amu=1.;
  CFreal Bmu=0.0165;
  CFreal Dmu=20.5;
  CFreal Aeps1=0.05;
  CFreal Aeps2=1.;
  CFreal Beps2=1.;
  CFreal Rt=K*K/(epsilon*m_MuLam);
  CFreal Ry=sqrt(K)*wallNearestDistance[cid]/m_MuLam;  
  CFreal fmu=(1.-Amu*exp(-Bmu*Ry))*(1.-Amu*exp(-Bmu*Ry))*(1.+Dmu/Rt);
  CFreal feps1=1.+(Aeps1/fmu)*(Aeps1/fmu)*(Aeps1/fmu);
  CFreal feps2=1.-Aeps2*exp(-(Beps2*Rt)*(Beps2*Rt));
*/

{
const std::vector<Node*>& nodes = cell.getNodes();
Node& node0 = *nodes[0];
Node& node1 = *nodes[1];
Node& node2 = *nodes[2];
Node& node3 = *nodes[3];
if (((node0[0]+node1[0]+node2[0]+node3[0])/4.>1.954)&&((node0[0]+node1[0]+node2[0]+node3[0])/4.<1.955)){
CFreal uplus=estate[1]/ustar;
CFreal upluslam=yplus;
CFreal upluslog=log(yplus)/0.41+5.1;
  CFLog(INFO,    wallNearestDistance[cid] << " "
              << wallNearestSegment[cid] << " " 
              << Rt << " " 
              << yplus << " " 
              << fmu << " " 
              << feps1 << " " 
              << feps2 << " " 
              << uplus << " " 
              << upluslam << " " 
              << upluslog << " \n" );   
}
}

  // calculate turbulent viscosity (muturb)
  const CFreal muturb      = fmu*RhoElm*m_Cmu*K*K/epsilon;
  const CFreal nuturb      = muturb/RhoElm;
  const CFreal nuk         = (m_MuLam+muturb/m_SigmaK)/RhoElm;
  const CFreal nueps       = (m_MuLam+muturb/m_SigmaEPS)/RhoElm;
  getMethodData().setMuElm(m_MuLam+muturb);

  CFuint Pi,Ui,Vi,Wi,Ki,EPSi,Pj,Uj,Vj,Wj,Kj,EPSj;
  CFreal val;
  CFreal S2 = 2.*((dudx*dudx)+(dvdy*dvdy)+(dwdz*dwdz))+((dudy+dvdx)*(dudy+dvdx)+(dudz+dwdx)*(dudz+dwdx)+(dvdz+dwdy)*(dvdz+dwdy));

  // SUPG
  CFreal umag=sqrt(estate[1]*estate[1]+estate[2]*estate[2]+estate[3]*estate[3]);
  CFreal tauSUPGK[3]={0.,0.,0.},tauSUPGEPS[3]={0.,0.,0.};
  const CellProps::CellData &celldata=m_cellprops.getCellData();

  if (umag!=0.)
  {
    CFreal sx=estate[1]/umag;
    CFreal sy=estate[2]/umag;
    CFreal sz=estate[3]/umag;
    CFreal h=0.;
    h+=fabs(sx*celldata.nx[0]+sy*celldata.ny[0]+sz*celldata.nz[0]);
    h+=fabs(sx*celldata.nx[1]+sy*celldata.ny[1]+sz*celldata.nz[1]);
    h+=fabs(sx*celldata.nx[2]+sy*celldata.ny[2]+sz*celldata.nz[2]);
    h+=fabs(sx*celldata.nx[3]+sy*celldata.ny[3]+sz*celldata.nz[3]);
    h/=(2.*celldata.vol);
    h=1./h;
    CFreal ree=umag*h/(2.*nuk);
    CFreal xi=max(0.,min(ree/3.,1.));
    CFreal tau=h*xi/(2.*umag);
    tauSUPGK[0]+=tau*estate[1];
    tauSUPGK[1]+=tau*estate[2];
    tauSUPGK[2]+=tau*estate[3];
    ree=umag*h/(2.*nueps);
    xi=max(0.,min(ree/3.,1.));
    tau=h*xi/(2.*umag);
    tauSUPGEPS[0]+=tau*estate[1];
    tauSUPGEPS[1]+=tau*estate[2];
    tauSUPGEPS[2]+=tau*estate[3];
  }

  // Nonlinear SUPG 
  CFreal uparK[3]={0.,0.,0.};
  CFreal normgradK[3]={0.,0.,0.};
  CFreal absgradK=sqrt(dKdx*dKdx+dKdy*dKdy+dKdz*dKdz);
  if (fabs(absgradK)>1.e-10){
    normgradK[0]=dKdx/absgradK;
    normgradK[1]=dKdy/absgradK;
    normgradK[2]=dKdz/absgradK;
    uparK[0]=(estate[1]*normgradK[0]+estate[2]*normgradK[1]+estate[3]*normgradK[2])*normgradK[0];
    uparK[1]=(estate[1]*normgradK[0]+estate[2]*normgradK[1]+estate[3]*normgradK[2])*normgradK[1];
    uparK[2]=(estate[1]*normgradK[0]+estate[2]*normgradK[1]+estate[3]*normgradK[2])*normgradK[2];
  } 
  CFreal uparKmag=sqrt(uparK[0]*uparK[0]+uparK[1]*uparK[1]+uparK[2]*uparK[2]);

  CFreal uparEPS[3]={0.,0.,0.};
  CFreal normgradEPS[3]={0.,0.,0.};
  CFreal absgradEPS=sqrt(dEPSdx*dEPSdx+dEPSdy*dEPSdy+dEPSdz*dEPSdz);
  if (fabs(absgradEPS)>1.e-10){
    normgradEPS[0]=dEPSdx/absgradEPS;
    normgradEPS[1]=dEPSdy/absgradEPS;
    normgradEPS[2]=dEPSdz/absgradEPS;
    uparEPS[0]=(estate[1]*normgradEPS[0]+estate[2]*normgradEPS[1]+estate[3]*normgradEPS[2])*normgradEPS[0];
    uparEPS[1]=(estate[1]*normgradEPS[0]+estate[2]*normgradEPS[1]+estate[3]*normgradEPS[2])*normgradEPS[1];
    uparEPS[2]=(estate[1]*normgradEPS[0]+estate[2]*normgradEPS[1]+estate[3]*normgradEPS[2])*normgradEPS[2];
  }
  CFreal uparEPSmag=sqrt(uparEPS[0]*uparEPS[0]+uparEPS[1]*uparEPS[1]+uparEPS[2]*uparEPS[2]);

  if (uparKmag!=0.)
  {
    CFreal sx=uparK[0]/uparKmag;
    CFreal sy=uparK[1]/uparKmag;
    CFreal sz=uparK[2]/uparKmag;
    CFreal h=0.;
    h+=fabs(sx*celldata.nx[0]+sy*celldata.ny[0]+sz*celldata.nz[0]);
    h+=fabs(sx*celldata.nx[1]+sy*celldata.ny[1]+sz*celldata.nz[1]);
    h+=fabs(sx*celldata.nx[2]+sy*celldata.ny[2]+sz*celldata.nz[2]);
    h+=fabs(sx*celldata.nx[3]+sy*celldata.ny[3]+sz*celldata.nz[3]);
    h/=(2.*celldata.vol);
    h=1./h;
    CFreal ree=uparKmag*h/(2.*nuk);
    CFreal xi=max(0.,min(ree/3.,1.));
    CFreal tau=h*xi/(2.*uparKmag);
    tauSUPGK[0]+=tau*uparK[0];
    tauSUPGK[1]+=tau*uparK[1];
    tauSUPGK[2]+=tau*uparK[2];
  }
  if (uparEPSmag!=0.)
  {
    CFreal sx=uparEPS[0]/uparEPSmag;
    CFreal sy=uparEPS[1]/uparEPSmag;
    CFreal sz=uparEPS[2]/uparEPSmag;
    CFreal h=0.;
    h+=fabs(sx*celldata.nx[0]+sy*celldata.ny[0]+sz*celldata.nz[0]);
    h+=fabs(sx*celldata.nx[1]+sy*celldata.ny[1]+sz*celldata.nz[1]);
    h+=fabs(sx*celldata.nx[2]+sy*celldata.ny[2]+sz*celldata.nz[2]);
    h+=fabs(sx*celldata.nx[3]+sy*celldata.ny[3]+sz*celldata.nz[3]);
    h/=(2.*celldata.vol);
    h=1./h;
    CFreal ree=uparEPSmag*h/(2.*nueps);
    CFreal xi=max(0.,min(ree/3.,1.));
    CFreal tau=h*xi/(2.*uparEPSmag);
    tauSUPGEPS[0]+=tau*uparEPS[0];
    tauSUPGEPS[1]+=tau*uparEPS[1];
    tauSUPGEPS[2]+=tau*uparEPS[2];
  }

//  tauSUPG*=0.;

  Pi=0; Ui=1; Vi=2; Wi=3; Ki=4; EPSi=5;
  for(CFuint i=0; i<3; ++i) {

    //CFreal u_ni   = estate[1]*celldata.nx[i]+estate[2]*celldata.ny[i];
    CFreal tauSUPGK_u_ni = tauSUPGK[0]*celldata.nx[i]+tauSUPGK[1]*celldata.ny[i]+tauSUPGK[2]*celldata.nz[i];
    CFreal tauSUPGEPS_u_ni = tauSUPGEPS[0]*celldata.nx[i]+tauSUPGEPS[1]*celldata.ny[i]+tauSUPGEPS[2]*celldata.nz[i];

    Pj=0; Uj=1; Vj=2; Wj=3; Kj=4; EPSj=5;
    for(CFuint j=0; j<3; ++j) {

      for(CFuint k=0; k<3; ++k){

        CFreal uk   = estate[1];
        CFreal vk   = estate[2];
        CFreal wk   = estate[3];
        CFreal uknj = uk*celldata.nx[j]+vk*celldata.ny[j]+wk*celldata.nz[j];

        // Convection (Standard + SUPG)
        val  = 1./60.*uknj*(1.+kronecker(i,k));
        val += tauSUPGK_u_ni/(36.*celldata.vol)*uknj;
        adata.A(Ki,Kj) += val;

        // Convection E (Standard + SUPG)
        val  = 1./60.*uknj*(1.+kronecker(i,k));
        val += tauSUPGEPS_u_ni/(36.*celldata.vol)*uknj;
        adata.A(EPSi,EPSj) += val;

      }

      //diffusion K (Standard)
      val = (nuk)/(9.*celldata.vol)*(celldata.nx[i]*celldata.nx[j]+celldata.ny[i]*celldata.ny[j]+celldata.nz[i]*celldata.nz[j]);
      adata.A(Ki,Kj) += val;

      //diffusion E (Standard)
      val = (nueps)/(9.*celldata.vol)*(celldata.nx[i]*celldata.nx[j]+celldata.ny[i]*celldata.ny[j]+celldata.nz[i]*celldata.nz[j]);
      adata.A(EPSi,EPSj) += val;

      // Time K (Standard)
      val  = celldata.vol/20.*(1.+kronecker(i,j));
      adata.T(Ki,Kj) += val;

      // Time E (Standard)
      val  = celldata.vol/20.*(1.+kronecker(i,j));
      adata.T(EPSi,EPSj) += val;

      // Production Term K
      val = -(nuturb*S2);
      val *=  celldata.vol/20.*(1.+kronecker(i,j));
      adata.b[Ki] += val;

      // Source Term K
      val = fmu*m_Cmu*K/nuturb;
      val *=  celldata.vol/20.*(1.+kronecker(i,j));
      adata.A(Ki,Kj) += val;

      // Source Term 1 E
      val = -feps1*m_Ceps1*fmu*m_Cmu*K*S2;
      val *= celldata.vol/20.*(1.+kronecker(i,j));
      adata.b[EPSi] += val;

      // Source Term 2 E
      val = feps2*m_Ceps2*m_Cmu*fmu*K/nuturb;
      val *= celldata.vol/20.*(1.+kronecker(i,j));
      adata.A(EPSi,EPSj) += val;
      
      // Chien Sources - D (into K)
      val = 2.*m_MuLam/(RhoElm*wallNearestDistance[cid]*wallNearestDistance[cid]);
      val *= celldata.vol/20.*(1.+kronecker(i,j));
      adata.A(Ki,Kj) += val;

      // Chien Sources - E (into EPS)
      val = 2.*m_MuLam/(RhoElm*wallNearestDistance[cid]*wallNearestDistance[cid])*exp(-0.5*yplus);
      val *= celldata.vol/20.*(1.+kronecker(i,j));
      adata.A(EPSi,EPSj) += val;

      Pj+=nbEqs; Uj+=nbEqs; Vj+=nbEqs; Wj+=nbEqs; Kj+=nbEqs; EPSj+=nbEqs;
    }
    Pi+=nbEqs; Ui+=nbEqs; Vi+=nbEqs; Wi+=nbEqs; Ki+=nbEqs; EPSi+=nbEqs;
  }

  return;

#undef kronecker

}

//////////////////////////////////////////////////////////////////////////////

void StandardKEpsilon::post ()
{
  CFAUTOTRACE;
  CFLog(INFO, getClassName() << "K Top limiting count:    " << m_KTopLimitCounter << "\n");
  CFLog(INFO, getClassName() << "K Bottom limiting count: " << m_KBottomLimitCounter << "\n");
  CFLog(INFO, getClassName() << "E Top limiting count:    " << m_ETopLimitCounter << "\n");
  CFLog(INFO, getClassName() << "E Bottom limiting count: " << m_EBottomLimitCounter << "\n");
}

//////////////////////////////////////////////////////////////////////////////

std::vector< Common::SafePtr< BaseDataSocketSource > > StandardKEpsilon::providesSockets()
{
  CFAUTOTRACE;
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr< BaseDataSocketSink > > StandardKEpsilon::needsSockets()
{
  CFAUTOTRACE;
  std::vector< SafePtr< BaseDataSocketSink > > result = UFEMTerm::needsSockets();

  result.push_back(&socket_interStates);
  result.push_back(&socket_wallNearestDistance);
  result.push_back(&socket_wallNearestSegment);
  result.push_back(&socket_wallNearestVelocityGradient);
  result.push_back(&socket_connBndFace2InnerCell);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    }  // namespace TetraP1P1Cell
  }  // namespace UFEM
}  // namespace COOLFluiD

