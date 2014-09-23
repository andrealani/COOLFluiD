#include "Environment/ObjectProvider.hh"

#include "UFEM/TetraP1P1Cell/SpalartAllmaras.hh"
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

Environment::ObjectProvider < SpalartAllmaras,
                              UFEMTerm,
                              UFEMTetraP1P1CellPlugin,
                              UFEMTerm::NARGS >
aTetraP1P1Cell_SpalartAllmaras_Provider ( "TetraP1P1Cell_SpalartAllmaras" );

//////////////////////////////////////////////////////////////////////////////

void SpalartAllmaras::defineConfigOptions ( Config::OptionList& options )
{
  CFAUTOTRACE;

  options.addConfigOption< CFreal >("Cb1",        "Cb1");
  options.addConfigOption< CFreal >("Cb2",        "Cb2");
  options.addConfigOption< CFreal >("Cw2",        "Cw2");
  options.addConfigOption< CFreal >("Cw3",        "Cw3");
  options.addConfigOption< CFreal >("Cv1",        "Cv1");
  options.addConfigOption< CFreal >("Cv2",        "Cv2");
  options.addConfigOption< CFreal >("K",          "Prandtl Constant");
  options.addConfigOption< CFreal >("MuLam",      "Laminar viscosity");
  options.addConfigOption< CFreal >("Sigma",      "Turbulent Prandtl Number");
  options.addConfigOption< CFreal >("TopLimit",   "Turb/Laminar viscosity ratio, above the turbulent viscosity is cropped.");
  options.addConfigOption< CFreal >("BottomLimit","Turb/Laminar viscosity ratio, below the turbulent viscosity is cropped.");
}


//////////////////////////////////////////////////////////////////////////////

SpalartAllmaras::SpalartAllmaras ( const std::string& name, Common::SafePtr<ElemProps> props  ) :
  UFEMTerm ( name, props ),
  m_cellprops ( *props.d_castTo<TetraP1P1Cell::CellProps>() ),
  socket_interStates("interStates"),
  socket_wallDistance("wallDistance")
{
  CFAUTOTRACE;
  addConfigOptionsTo(this);

  //Constants Default Values
  m_Cb1         = 0.1355;
  m_Cb2         = 0.622;
  m_Cw2         = 0.3;
  m_Cw3         = 2.;
  m_Cv1         = 7.1;
  m_Cv2         = 5.;
  m_K           = 0.41;
  m_Sigma       = 2./3.;
  m_MuLam       = 1.;
  m_TopLimit    = 1.e+5;
  m_BottomLimit = 1.e-5;

  setParameter( "Cb1",          &m_Cb1 );
  setParameter( "Cb2",          &m_Cb2 );
  setParameter( "Cw2",          &m_Cw2 );
  setParameter( "Cw3",          &m_Cw3 );
  setParameter( "Cv1",          &m_Cv1 );
  setParameter( "Cv2",          &m_Cv2 );
  setParameter( "K",            &m_K );
  setParameter( "Sigma",        &m_Sigma );
  setParameter( "MuLam",        &m_MuLam );
  setParameter( "TopLimit",     &m_TopLimit );
  setParameter( "BottomLimit",  &m_BottomLimit );
}

//////////////////////////////////////////////////////////////////////////////

SpalartAllmaras::~SpalartAllmaras() 
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void SpalartAllmaras::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;
  UFEMTerm::configure ( args );

  CFLog(INFO, getClassName() << ": Cb1: "                                      << m_Cb1         << "\n");
  CFLog(INFO, getClassName() << ": Cb2: "                                      << m_Cb2         << "\n");
  CFLog(INFO, getClassName() << ": Cw2: "                                      << m_Cw2         << "\n");
  CFLog(INFO, getClassName() << ": Cw3: "                                      << m_Cw3         << "\n");
  CFLog(INFO, getClassName() << ": Cv1: "                                      << m_Cv1         << "\n");
  CFLog(INFO, getClassName() << ": Cv2: "                                      << m_Cv2         << "\n");
  CFLog(INFO, getClassName() << ": K (Prandtl Constant): "                     << m_K           << "\n");
  CFLog(INFO, getClassName() << ": Sigma (Turbulent Prandtl Number): "         << m_Sigma       << "\n");
  CFLog(INFO, getClassName() << ": MuLam (Laminar viscosity): "                << m_MuLam       << "\n");
  CFLog(INFO, getClassName() << ": TopLimit (viscosity limit by ratio): "      << m_TopLimit    << "\n");
  CFLog(INFO, getClassName() << ": BottomLimit (viscosity limit by ratio): "   << m_BottomLimit << "\n");
}

//////////////////////////////////////////////////////////////////////////////

void SpalartAllmaras::setup ()
{
  CFAUTOTRACE;
  UFEMTerm::setup ();

  socket_interStates.setParentNamespace ( getMethodData().getNamespace() );
  socket_wallDistance.setParentNamespace ( getMethodData().getNamespace() );
}

//////////////////////////////////////////////////////////////////////////////

void SpalartAllmaras::unsetup ()
{
  CFAUTOTRACE;
  UFEMTerm::unsetup ();
}

//////////////////////////////////////////////////////////////////////////////

void SpalartAllmaras::pre ()
{
  CFAUTOTRACE;
  m_BottomLimitCounter=0;
  m_TopLimitCounter=0;
}

//////////////////////////////////////////////////////////////////////////////

void SpalartAllmaras::compute ( const Framework::GeometricEntity& cell , AssemblyData& adata )
{
  CFAUTOTRACE;


#define kronecker(i,j) ((i)==(j)?1.:0.)

  CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CellProps::CellData &celldata=m_cellprops.getCellData();

  RealVector estate(nbEqs);
  estate = 0.;

  const vector< State* >& states = cell.getStates();
  DataHandle<Framework::State*> interStates = socket_interStates.getDataHandle();
  DataHandle<CFreal> wallDistance = socket_wallDistance.getDataHandle();

  CFreal dudx=0.,dudy=0.,dudz=0.,dvdx=0.,dvdy=0.,dvdz=0.,dwdx=0.,dwdy=0.,dwdz=0.,dnuturbdx=0.,dnuturbdy=0.,dnuturbdz=0.,WD=0.;

  const CFreal invvol=1./(3.*m_cellprops.getCellData().vol);
  const CFreal* nx=celldata.nx;
  const CFreal* ny=celldata.ny;
  const CFreal* nz=celldata.nz;

  const CFreal toplimit = (m_MuLam/getMethodData().getRhoElm())*m_TopLimit;
  const CFreal bottomlimit = (m_MuLam/getMethodData().getRhoElm())*m_BottomLimit;

  for (CFuint iState=0; iState<4; ++iState) {
    State& state = *states[iState];
    State& interState = *interStates[state.getLocalID()];
    const CFuint nodeID = state.getLocalID();
    const CFreal walld = wallDistance[nodeID];
    if ((interState[4]<bottomlimit)&&(walld>1e-15)) {
      interState[4]=bottomlimit;
      //state[4]=bottomlimit;
      //cout << "M " << flush;
      m_BottomLimitCounter++;
    }
    if ((interState[4]>toplimit)&&(walld>1e-15)) {
      interState[4]=toplimit;
      //state[4]=toplimit;
      //cout << "P " << flush;
      m_TopLimitCounter++;
    }

    dudx            +=interState[1]*nx[iState];
    dudy            +=interState[1]*ny[iState];
    dudz            +=interState[1]*nz[iState];
    dvdx            +=interState[2]*nx[iState];
    dvdy            +=interState[2]*ny[iState];
    dvdz            +=interState[2]*nz[iState];
    dwdx            +=interState[3]*nx[iState];
    dwdy            +=interState[3]*ny[iState];
    dwdz            +=interState[3]*nz[iState];
    dnuturbdx       +=interState[4]*nx[iState];
    dnuturbdy       +=interState[4]*ny[iState];
    dnuturbdz       +=interState[4]*nz[iState];
    WD              +=walld;
    for (CFuint iEq=0; iEq<nbEqs; ++iEq) {
      estate[iEq] += interState[iEq];
    }
  }
  dudx         *= invvol;
  dudy         *= invvol;
  dudz         *= invvol;
  dvdx         *= invvol;
  dvdy         *= invvol;
  dvdz         *= invvol;
  dwdx         *= invvol;
  dwdy         *= invvol;
  dwdz         *= invvol;
  dnuturbdx    *= invvol;
  dnuturbdy    *= invvol;
  dnuturbdz    *= invvol;
  for (CFuint iEq=0; iEq<nbEqs; ++iEq) estate[iEq] *= 1./4.;
  WD *= 1./3.;
  CFreal gradnuturb[3]={0.,0.,0.};
         gradnuturb[0]=dnuturbdx;
         gradnuturb[1]=dnuturbdy;
         gradnuturb[2]=dnuturbdz;

////////////////////////////////////

  // General physichal quantities
  const CFreal RhoElm = getMethodData().getRhoElm();
  const CFreal nulam  = m_MuLam/RhoElm;
  CFreal nuturb = estate[4];
  const CFreal muturb = nuturb*RhoElm;
  const CFreal S      = sqrt((dudy-dvdx)*(dudy-dvdx)+(dudz-dwdx)*(dudz-dwdx)+(dvdz-dwdy)*(dvdz-dwdy));
  const CFreal chi    = muturb/m_MuLam;
  const CFreal fv1    = (chi*chi*chi)/(chi*chi*chi+m_Cv1*m_Cv1*m_Cv1);
//  // original SA
//  const CFreal fv2    = 1.-chi/(1.+chi*fv1);
//  const CFreal Stilde = S+nuturb/(m_K*m_K*WD*WD)*fv2;
  // modification proposed by Spalart to avoid negative Stilde
  const CFreal invfv2 = 1.+chi/m_Cv2;
  const CFreal fv2    = 1./(invfv2*invfv2*invfv2);
  const CFreal fv3    = (1.+chi*fv1)*(1.-fv2)/chi;
  const CFreal Stilde = S*fv3+nuturb/(m_K*m_K*WD*WD)*fv2;
//  const CFreal r      = nuturb/(Stilde*m_K*m_K*WD*WD);
  const CFreal r      = min(nuturb/(Stilde*m_K*m_K*WD*WD),10.);
  const CFreal g      = r+m_Cw2*(r*r*r*r*r*r-r);
  const CFreal cw36   = m_Cw3*m_Cw3*m_Cw3*m_Cw3*m_Cw3*m_Cw3;
//  const CFreal fw     = g*pow((1.+cw36)/(g*g*g*g*g*g+cw36),1./6.);
  const CFreal fw     = pow((1.+cw36)/(1.+cw36/(g*g*g*g*g*g)),1./6.);
  const CFreal Cw1    = m_Cb1/(m_K*m_K)+((1.+m_Cb2)/m_Sigma);

  // Set Element viscosity (Laminar + Turbulent)
  getMethodData().setMuElm(m_MuLam+fv1*muturb);

  // SUPG
  CFreal umag=sqrt(estate[1]*estate[1]+estate[2]*estate[2]+estate[3]*estate[3]);
  CFreal tauSUPG[3]={0.,0.,0.};
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
    CFreal ree=umag*h/(2.*nuturb);
    CFreal xi=max(0.,min(ree/3.,1.));
    CFreal tau=h*xi/(2.*umag);
    tauSUPG[0]+=tau*estate[1];
    tauSUPG[1]+=tau*estate[2];
    tauSUPG[2]+=tau*estate[3];

  }

  // Nonlinear SUPG 
  CFreal upar[3]={0.,0.,0.};
  CFreal normgradnuturb[3]={0.,0.,0.};
  CFreal absgradnuturb=sqrt(dnuturbdx*dnuturbdx+dnuturbdy*dnuturbdy+dnuturbdz*dnuturbdz);
  if (fabs(absgradnuturb)>1.e-10){
    normgradnuturb[0]=dnuturbdx/absgradnuturb;
    normgradnuturb[1]=dnuturbdy/absgradnuturb;
    normgradnuturb[2]=dnuturbdz/absgradnuturb;
    upar[0]=(estate[1]*normgradnuturb[0]+estate[2]*normgradnuturb[1]+estate[3]*normgradnuturb[2])*normgradnuturb[0];
    upar[1]=(estate[1]*normgradnuturb[0]+estate[2]*normgradnuturb[1]+estate[3]*normgradnuturb[2])*normgradnuturb[1];
    upar[2]=(estate[1]*normgradnuturb[0]+estate[2]*normgradnuturb[1]+estate[3]*normgradnuturb[2])*normgradnuturb[2];
  } 
  CFreal uparmag=sqrt(upar[0]*upar[0]+upar[1]*upar[1]+upar[2]*upar[2]);

  if (uparmag!=0.)
  {
    CFreal sx=upar[0]/uparmag;
    CFreal sy=upar[1]/uparmag;
    CFreal sz=upar[2]/uparmag;
    CFreal h=0.;
    h+=fabs(sx*celldata.nx[0]+sy*celldata.ny[0]+sz*celldata.nz[0]);
    h+=fabs(sx*celldata.nx[1]+sy*celldata.ny[1]+sz*celldata.nz[1]);
    h+=fabs(sx*celldata.nx[2]+sy*celldata.ny[2]+sz*celldata.nz[2]);
    h+=fabs(sx*celldata.nx[3]+sy*celldata.ny[3]+sz*celldata.nz[3]);
    h/=(2.*celldata.vol);
    h=1./h;
    CFreal ree=uparmag*h/(2.*nuturb);
    CFreal xi=max(0.,min(ree/3.,1.));
    CFreal tau=h*xi/(2.*uparmag);
    tauSUPG[0]+=tau*upar[0];
    tauSUPG[1]+=tau*upar[1];
    tauSUPG[2]+=tau*upar[2];
  }

//  tauSUPG*=0.;

  CFuint Pi,Ui,Vi,Wi,NUTURBi,Pj,Uj,Vj,Wj,NUTURBj;
  CFreal val;

  Pi=0; Ui=1; Vi=2; Wi=3; NUTURBi=4;
  for(CFuint i=0; i<3; ++i) {

    //CFreal u_ni   = estate[1]*celldata.nx[i]+estate[2]*celldata.ny[i];
    CFreal tauSUPG_u_ni = tauSUPG[0]*celldata.nx[i]+tauSUPG[1]*celldata.ny[i];

    Pj=0; Uj=1; Vj=2; Wj=3; NUTURBj=4;
    for(CFuint j=0; j<3; ++j) {

      for(CFuint k=0; k<3; ++k){

        CFreal uk   = estate[1];
        CFreal vk   = estate[2];
        CFreal wk   = estate[3];
        CFreal uknj = uk*celldata.nx[j]+vk*celldata.ny[j]+wk*celldata.nz[j];

        // Convection (Standard + SUPG)
        val  = 1./60.*uknj*(1.+kronecker(i,k));
        val += tauSUPG_u_ni/(36.*celldata.vol)*uknj;
        adata.A(NUTURBi,NUTURBj) += val;

      }

      // Diffusion (Standard)
      val = ((nulam+nuturb)/m_Sigma);
      val *= 1./(9.*celldata.vol)*(celldata.nx[i]*celldata.nx[j]+celldata.ny[i]*celldata.ny[j]+celldata.nz[i]*celldata.nz[j]);
      adata.A(NUTURBi,NUTURBj) += val;


      // Time (Standard)
      val  = celldata.vol/20.*(1.+kronecker(i,j));
      adata.T(NUTURBi,NUTURBj) += val;

      // Addition to diffusion term (Standard)
      val = -(m_Cb2/m_Sigma)*(gradnuturb[0]*gradnuturb[0]+gradnuturb[1]*gradnuturb[1])*celldata.vol/20.*(1.+kronecker(i,j));
      adata.b[NUTURBi] += val;

      // Production term (Standard)
      val = -(m_Cb1*Stilde)*celldata.vol/20.*(1.+kronecker(i,j));
      adata.A(NUTURBi,NUTURBj) += val;

      // Destruction Term (Standard)
      val = fw*Cw1*nuturb/(WD*WD)*celldata.vol/20.*(1.+kronecker(i,j));
      adata.A(NUTURBi,NUTURBj) += val;

      Pj+=nbEqs; Uj+=nbEqs; Vj+=nbEqs; Wj+=nbEqs; NUTURBj+=nbEqs;
    }
    Pi+=nbEqs; Ui+=nbEqs; Vi+=nbEqs; Wi+=nbEqs; NUTURBi+=nbEqs;
  }

  return;

#undef kronecker

}

//////////////////////////////////////////////////////////////////////////////

void SpalartAllmaras::post ()
{
  CFAUTOTRACE;
  CFLog(INFO, getClassName() << "Top limiting count:    " << m_TopLimitCounter << "\n");
  CFLog(INFO, getClassName() << "Bottom limiting count: " << m_BottomLimitCounter << "\n");
}

//////////////////////////////////////////////////////////////////////////////

std::vector< SafePtr< BaseDataSocketSink > > SpalartAllmaras::needsSockets()
{
  CFAUTOTRACE;
  std::vector< SafePtr< BaseDataSocketSink > > result = UFEMTerm::needsSockets();

  result.push_back(&socket_interStates);
  result.push_back(&socket_wallDistance);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    }  // namespace TetraP1P1Cell
  }  // namespace UFEM
}  // namespace COOLFluiD

