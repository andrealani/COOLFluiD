#include "Environment/ObjectProvider.hh"

#include "UFEM/TriagP1P1Cell/SpalartAllmaras.hh"
#include "UFEM/TriagP1P1Cell/UFEMTriagP1P1Cell.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {
  namespace UFEM {
    namespace TriagP1P1Cell {

using namespace std;
using namespace COOLFluiD::Environment;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider < SpalartAllmaras,
                              UFEMTerm,
                              UFEMTriagP1P1CellPlugin,
                              UFEMTerm::NARGS >
aTriagP1P1Cell_SpalartAllmaras_Provider ( "TriagP1P1Cell_SpalartAllmaras" );

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
  m_cellprops ( *props.d_castTo<TriagP1P1Cell::CellProps>() ),
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

  CFreal dudx=0.,dudy=0.,dvdx=0.,dvdy=0.,dnuturbdx=0.,dnuturbdy=0.,WD=0.;

  const CFreal invvol=1./(2.*m_cellprops.getCellData().vol);
  const CFreal* nx=celldata.nx;
  const CFreal* ny=celldata.ny;

  const CFreal toplimit = (m_MuLam/getMethodData().getRhoElm())*m_TopLimit;
  const CFreal bottomlimit = (m_MuLam/getMethodData().getRhoElm())*m_BottomLimit;

  for (CFuint iState=0; iState<3; ++iState) {
    State& state = *states[iState];
    State& interState = *interStates[state.getLocalID()];
    const CFuint nodeID = state.getLocalID();
    const CFreal walld = wallDistance[nodeID];
    if ((interState[3]<bottomlimit)&&(walld>1e-15)) {
      interState[3]=bottomlimit;
      //state[3]=bottomlimit;
      //cout << "M " << flush;
      m_BottomLimitCounter++;
    }
    if ((interState[3]>toplimit)&&(walld>1e-15)) {
      interState[3]=toplimit;
      //state[3]=toplimit;
      //cout << "P " << flush;
      m_TopLimitCounter++;
    }

    dudx            +=interState[1]*nx[iState];
    dudy            +=interState[1]*ny[iState];
    dvdx            +=interState[2]*nx[iState];
    dvdy            +=interState[2]*ny[iState];
    dnuturbdx       +=interState[3]*nx[iState];
    dnuturbdy       +=interState[3]*ny[iState];
    WD              +=walld;
    for (CFuint iEq=0; iEq<nbEqs; ++iEq) {
      estate[iEq] += interState[iEq];
    }
  }
  dudx         *= invvol;
  dudy         *= invvol;
  dvdx         *= invvol;
  dvdy         *= invvol;
  dnuturbdx    *= invvol;
  dnuturbdy    *= invvol;
  for (CFuint iEq=0; iEq<nbEqs; ++iEq) estate[iEq] *= 1./3.;
  WD *= 1./3.;
  CFreal gradnuturb[2]={0.,0.};
         gradnuturb[0]=dnuturbdx;
         gradnuturb[1]=dnuturbdy;

////////////////////////////////////

  // General physichal quantities
  const CFreal RhoElm = getMethodData().getRhoElm();
  const CFreal nulam  = m_MuLam/RhoElm;
  const CFreal nuturb = estate[3];
  const CFreal muturb = nuturb*RhoElm;
  const CFreal S      = sqrt((dudy-dvdx)*(dudy-dvdx));
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
  CFreal umag=sqrt(estate[1]*estate[1]+estate[2]*estate[2]);
  CFreal tauSUPG[2]={0.,0.};
  if (umag!=0.)
  {
    CFreal sx=estate[1]/umag;
    CFreal sy=estate[2]/umag;
    CFreal h=0.;
    h+=fabs(sx*celldata.nx[0]+sy*celldata.ny[0]);
    h+=fabs(sx*celldata.nx[1]+sy*celldata.ny[1]);
    h+=fabs(sx*celldata.nx[2]+sy*celldata.ny[2]);
    h/=(2.*celldata.vol);
    h=1./h;
    CFreal ree=umag*h/(2.*nuturb);
    CFreal xi=max(0.,min(ree/3.,1.));
    CFreal tau=h*xi/(2.*umag);
    tauSUPG[0]+=tau*estate[1];
    tauSUPG[1]+=tau*estate[2];

  }

  // Nonlinear SUPG 
  CFreal upar[2]={0.,0.};
  CFreal normgradnuturb[2]={0.,0.};
  CFreal absgradnuturb=sqrt(dnuturbdx*dnuturbdx+dnuturbdy*dnuturbdy);
  if (fabs(absgradnuturb)>1.e-10){
    normgradnuturb[0]=dnuturbdx/absgradnuturb;
    normgradnuturb[1]=dnuturbdy/absgradnuturb;
    upar[0]=(estate[1]*normgradnuturb[0]+estate[2]*normgradnuturb[1])*normgradnuturb[0];
    upar[1]=(estate[1]*normgradnuturb[0]+estate[2]*normgradnuturb[1])*normgradnuturb[1];
  } 
  CFreal uparmag=sqrt(upar[0]*upar[0]+upar[1]*upar[1]);
  if (uparmag!=0.)
  {
    CFreal sx=upar[0]/uparmag;
    CFreal sy=upar[1]/uparmag;
    CFreal h=0.;
    h+=fabs(sx*celldata.nx[0]+sy*celldata.ny[0]);
    h+=fabs(sx*celldata.nx[1]+sy*celldata.ny[1]);
    h+=fabs(sx*celldata.nx[2]+sy*celldata.ny[2]);
    h/=(2.*celldata.vol);
    h=1./h;
    CFreal ree=uparmag*h/(2.*nuturb);
    CFreal xi=max(0.,min(ree/3.,1.));
    CFreal tau=h*xi/(2.*uparmag);
    tauSUPG[0]+=tau*upar[0];
    tauSUPG[1]+=tau*upar[1];
  }

//  tauSUPG*=0.;

  CFuint Pi,Ui,Vi,NUTURBi,Pj,Uj,Vj,NUTURBj;
  CFreal val;

  Pi=0; Ui=1; Vi=2; NUTURBi=3;
  for(CFuint i=0; i<3; ++i) {

    //CFreal u_ni   = estate[1]*celldata.nx[i]+estate[2]*celldata.ny[i];
    CFreal tauSUPG_u_ni = tauSUPG[0]*celldata.nx[i]+tauSUPG[1]*celldata.ny[i];

    Pj=0; Uj=1; Vj=2; NUTURBj=3;
    for(CFuint j=0; j<3; ++j) {

      for(CFuint k=0; k<3; ++k){

        CFreal uk   = estate[1];
        CFreal vk   = estate[2];
        CFreal uknj = uk*celldata.nx[j]+vk*celldata.ny[j];

        // Convection (Standard + SUPG)
        val  = 1./24.*uknj*(1.+kronecker(i,k));
        val += tauSUPG_u_ni/(12.*celldata.vol)*uknj;
        adata.A(NUTURBi,NUTURBj) += val;

      }

      // Diffusion (Standard)
      val = ((nulam+nuturb)/m_Sigma)/(4.*celldata.vol)*(celldata.nx[i]*celldata.nx[j]+celldata.ny[i]*celldata.ny[j]);
      adata.A(NUTURBi,NUTURBj) += val;


      // Time (Standard)
      val  = celldata.vol/12.*(1.+kronecker(i,j));
      adata.T(NUTURBi,NUTURBj) += val;

      // Addition to diffusion term (Standard)
      val = -(m_Cb2/m_Sigma)*(gradnuturb[0]*gradnuturb[0]+gradnuturb[1]*gradnuturb[1])*celldata.vol/12.*(1.+kronecker(i,j));
      adata.b[NUTURBi] += val;

      // Production term (Standard)
      val = -(m_Cb1*Stilde)*celldata.vol/12.*(1.+kronecker(i,j));
      adata.A(NUTURBi,NUTURBj) += val;

      // Destruction Term (Standard)
      val = fw*Cw1*nuturb/(WD*WD)*celldata.vol/12.*(1.+kronecker(i,j));
      adata.A(NUTURBi,NUTURBj) += val;

      Pj+=nbEqs; Uj+=nbEqs; Vj+=nbEqs; NUTURBj+=nbEqs;
    }
    Pi+=nbEqs; Ui+=nbEqs; Vi+=nbEqs; NUTURBi+=nbEqs;
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

    }  // namespace TriagP1P1Cell
  }  // namespace UFEM
}  // namespace COOLFluiD

