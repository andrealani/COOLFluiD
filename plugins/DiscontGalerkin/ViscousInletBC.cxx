#include "Framework/MethodCommandProvider.hh"
#include "Framework/CFSide.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"
#include "Framework/FaceToCellGEBuilder.hh"

#include "DiscontGalerkin/ViscousInletBC.hh"
#include "DiscontGalerkin/DiscontGalerkin.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace DiscontGalerkin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ViscousInletBC,DiscontGalerkinSolverData,DiscontGalerkinModule >
  viscousInletBCProvider("ViscousInletBC");

//////////////////////////////////////////////////////////////////////////////

void ViscousInletBC::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >("Vars","Definition of the Variables.");
  options.addConfigOption< std::vector<std::string> >("Def","Definition of the Functions.");
}

//////////////////////////////////////////////////////////////////////////////

ViscousInletBC::ViscousInletBC(const std::string& name)
: ViscousBaseSolve(name),
  socket_rhs("rhs"),
  kappa1(0.4)

{
  addConfigOptionsTo(this);
  _functions = std::vector<std::string>();
  setParameter("Def",&_functions);

  _vars = std::vector<std::string>();
  setParameter("Vars",&_vars);
}

//////////////////////////////////////////////////////////////////////////////

ViscousInletBC::~ViscousInletBC()
{
}

//////////////////////////////////////////////////////////////////////////////

void ViscousInletBC::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  DiscontGalerkinSolverCom::configure(args);

  _vFunction.setFunctions(_functions);
  _vFunction.setVariables(_vars);
  try {
    _vFunction.parse();
  }
  catch (Common::ParserException& e) {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }
}

//////////////////////////////////////////////////////////////////////////////

void ViscousInletBC::setup()
{
  CFAUTOTRACE;
  ViscousBaseSolve::setup();
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  m_theta = getMethodData().getTheta();
  m_sigma = getMethodData().getSigma();

  m_kMatrix.resize(nbDim);
  for (CFuint k = 0; k < nbDim ; k++)
  {
    m_kMatrix[k].resize(nbDim);
    for (CFuint s = 0; s < nbDim ; s++)
    {
      m_kMatrix[k][s].resize(nbEqs,nbEqs);
      m_kMatrix[k][s]=0.;
    }
  }
  // set a pointer to the cells
  m_cells.reset(MeshDataStack::getActive()->getTrs("InnerCells"));

  SafePtr<vector<ElementTypeData> > elementType = MeshDataStack::getActive()->getElementTypeData();

  // rhs storage
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  const CFuint nbElemTypes = elementType->size();

  // loop over types since it can happen to deal with an hybrid mesh
  for (CFuint iType = 0; iType < nbElemTypes; ++iType) {
    const CFuint nbStatesInType = (*elementType)[iType].getNbStates();

      BlockAccumulator* ptr = getMethodData().getLinearSystemSolver()[0]->
       createBlockAccumulator(nbStatesInType,nbStatesInType,nbEqs);
      RealVector* vec = new RealVector(nbStatesInType*nbEqs);
      RealVector* vec2 = new RealVector(nbEqs);
      RealMatrix* mat = new RealMatrix(nbStatesInType*nbEqs,nbStatesInType*nbEqs);
      vector<RealVector>* residual = new vector<RealVector>(nbStatesInType,*vec2);

      DGElemTypeData elemTypeData(ptr,mat,vec,residual);

      m_mapElemData.insert(nbStatesInType, elemTypeData);
  }

  m_mapElemData.sortKeys();

}

//////////////////////////////////////////////////////////////////////////////


void ViscousInletBC::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

void ViscousInletBC::ExactRiemannSolver2D(Framework::State UL, Framework::State UR, RealVector S, Framework::State *US)
{
  CFAUTOTRACE;
// CFout << " left " << UL[0] << " PR " << UL[1] << " CL " << UL[2] << " CR " << UL[3] << CFendl; 
// CFout << " right " << UR[0] << " PR " << UR[1] << " CL " << UR[2] << " CR " << UR[3] << CFendl; 
//          Calculates the 2D Euler flux using the exact Riemann Solver
//          -----------------------------------------------------------

//          NB:  Input parameters are the left (UL) and right (UR) state
//          vectors of conserved variables,  with
//          U = (rho, rho u, rho v, rho w, rho E), and the side-area
//          vector, S, where S = d(s1,s2), where d is
//          the area of the cell face and (s1,s2) is the outward
//          pointing cell-face normal vector, pointing from cell L
//          to cell R.
//          The solution is rotated to the S vector frame of reference.
//          The EXACTRS subroutine also needs the xphi and PRESSURE subroutines
//          listed below.
//          The flux is returned in vector US.

    CFint I;
    CFreal SM,ut,vt,uu,vv,sv,pi,qi,ri,RSTAR,CSTAR,EEXP,SEXP,ei;

    CFreal RGAMMA=1.4;
    CFuint conv = 0;
    CFreal eps1 = 1.0E-6;
    CFreal eps2 = 1.0E-6;
    CFreal alpha = 1.0E0;
    CFreal d = std::sqrt((S[1]*S[1]+S[0]*S[0]));
    CFreal rdd = 1.0E0/(d+1.0E-20);
    CFreal s1 = S[0]*rdd;
    CFreal s2 = S[1]*rdd;
    CFreal Rrl = 1.0E0/UL[0];
    CFreal Rrr = 1.0E0/UR[0];
    CFreal ql = (s1*UL[1]+s2*UL[2])*Rrl;
    CFreal qr = (s1*UR[1]+s2*UR[2])*Rrr;

//     print*,'##',s1,s2,Rrl,Rrr,ql,qr;

   CFreal pl = pressure(UL);
   CFreal pr = pressure(UR);
   CFreal cl = std::sqrt((RGAMMA*pl*Rrl));
   CFreal cr = std::sqrt((RGAMMA*pr*Rrr));
   CFreal ps = (pl+pr)*0.5;

// CFout << " PL " << pl << " PR " << UL[1] << " CL " << UL[2] << " CR " << UL[3] << CFendl; 

//     print*,'##',pl,pr,cl,cr,ps,ql;

//     !    START OF EXACT RIEMANN SOLVER
   CFreal rconst = UR[0]*std::sqrt(pr);
   CFreal lconst = UL[0]*std::sqrt(pl);
   CFreal MR = rconst*xphi(ps/pr);
   CFreal ML = lconst*xphi(ps/pl);

//     !     Godunov's iteration
  while (conv == 0)
  {
    I=0;
    while ((I <20) && (conv == 0))
    {
      CFreal pn = (pr/MR+pl/ML+ql-qr)/(1.0E0/ML+1.0E0/MR);
      pn = alpha*std::max(eps1,pn)+(1.0E0-alpha)*ps;
      CFreal newmr = rconst*xphi(pn/pr);
      CFreal newml = lconst*xphi(pn/pl);
      CFreal t1 = CFreal(std::abs(newmr-MR));
      CFreal t2 = CFreal(std::abs(newml-ML));
      if ((max(t1,t2) < eps1) || (alpha < eps2)) conv=1;
      ps = pn;
      MR = newmr;
      ML = newml;
      I = I+1;
    }
    alpha = alpha*0.8;
  }
  SM = (pl-pr+ML*ql+MR*qr)/(ML+MR);

//     !print*,'##@!!',SM,pl,pr,ML,ql,MR,qr

//     !    Now find what is happening at the interface
  if (SM < 0.0)
  {
//        !      RIGHT OF THE CONTACT
    ut = UR[1]*Rrr-qr*s1;
    vt = UR[2]*Rrr-qr*s2;
    if ((ps/pr) > 1.0E0)
    {
//           !      RIGHT TRAVELLING SHOCK
      sv = qr+MR*Rrr;
      if (sv < 0.0E0)
      {
//              !         SUPERSONIC FROM RIGHT TO LEFT
        ri = UR[0];
        qi = qr;
        pi = pr;
      }
      else
      {
//              !         BETWEEN CONTACT AND RIGHT MOVING SHOCK
        ri = MR/(sv-SM);
        qi = SM;
        pi = ps;
      }
    }
    else
    {
//           !      RIGHT TRAVELLING EXPANSION
      EEXP=(qr+cr);
      if (EEXP < 0.0)
      {
//              !         SUPERSONIC FROM RIGHT TO LEFT
        ri = UR[0];
        qi = qr;
        pi = pr;
      }
      else
      {
        CFreal RSTAR=pow((ps*(std::pow(UR[0],RGAMMA))/pr),(1.0E0/RGAMMA));
        CFreal CSTAR=RGAMMA*ps/sqrt(RSTAR);
        CFreal SEXP=(CSTAR+SM);
        if (SEXP > 0.0E0)
        {
//                 !            BETWEEN CONTACT AND START OF EXPANSION
          ri = RSTAR;
          qi = SM;
          pi = ps;
        }
        else
        {
//                 !            IN THE EXPANSION WAVE
          qi=(((RGAMMA-1.0E0)*qr*0.5)-cr)* 2.0E0/(RGAMMA+1.0E0);
          ri=1.0E0+(RGAMMA-1.0E0)*(qi-qr)/(2.0E0*cr);
          ri=UR[0]*ri*(2.0E0/(RGAMMA-1.0E0));
          pi=(pr/(pow(UR[0],RGAMMA)))*(pow(ri,RGAMMA));
        }
      }
    }
  }
  else
  {
//        !      LEFT OF THE CONTACT
    ut = UL[1]*Rrl-ql*s1;
    vt = UL[2]*Rrl-ql*s2;
    if (ps/pl > 1.0E0)
    {
//           !        LEFT TRAVELLING SHOCK
      sv = ql-ML*Rrl;
      if (sv > 0.0E0)
      {
//              !          SUPERSONIC FROM LEFT TO RIGHT
        qi = ql;
        ri = UL[0];
        pi = pl;
      }
      else
      {
//              !          BETWEEN SHOCK AND CONTACT
        pi = ps;
        qi = SM;
        ri =  ML/(SM-sv);
      }
    }
    else
    {
//           !       LEFT TRAVELLING EXPANSION
      SEXP=ql-cl;
      if (SEXP > 0.0E0)
      {
//              !          Supersonic from left to right
        ri = UL[0];
        pi = pl;
        qi = ql;
      }
      else
      {
        RSTAR=std::pow(ps*( std::pow(UL[0] , RGAMMA) )/pl , 1.0/RGAMMA);
        CSTAR=RGAMMA*ps/std::sqrt(RSTAR);
        EEXP=SM-CSTAR;
//             !       if(ielem .eq. 79) write (*,*) 'SEXP=',SEXP,EEXP
        if (EEXP > 0.0E0)
        {
//                 !             In the expansion wave
          qi=((RGAMMA-1.0E0)*ql*0.5+cl)*2.0E0/(RGAMMA+1.0E0);
          ri=1.0E0-((RGAMMA-1.0E0)/(2.0E0*cl))*(qi-ql);
          ri=UL[0]*std::pow(ri,(2.0E0/(RGAMMA-1.0E0)));
          pi=(pl/(std::pow(UL[0],RGAMMA)))*(pow(ri,RGAMMA));
        }
        else
        {
//                 !             BETWEEN EXPANSION AND CONTACT
          ri = RSTAR;
          qi = SM;
          pi = ps;
//           !    if(ielem .eq. 79) write (*,*) 'ri,qi,pi=',ri,qi,pi
        }
      }
    }
  }


  uu = s1*qi+ut;
  vv = s2*qi+vt;

//     !     if(ielem .eq. 79) write (*,*) ri,uu,vv,pi,qi,'&&&&&'
//     !     Find the flux at the interface and scale by the face area, d
  (*US)[0] = ri;
  (*US)[1] = ri*uu;
  (*US)[2] = ri*vv;
  ei = pi/(RGAMMA-1.0E0)+0.5*ri*(uu*uu+vv*vv);
  (*US)[3] = ei;

//  CFout << " b0 " << US[0]<< " b1 " << US[1]<< " b2 " << US[2]<< " b3 " << US[3] << CFendl;
}
//////////////////////////////////////////////////////////////////////////////

void ViscousInletBC::ExactRiemannSolver3D(Framework::State UL, Framework::State UR, RealVector S, Framework::State *US)
{
  CFAUTOTRACE;
// CFout << " left " << UL[0] << " PR " << UL[1] << " CL " << UL[2] << " CR " << UL[3] << CFendl; 
// CFout << " right " << UR[0] << " PR " << UR[1] << " CL " << UR[2] << " CR " << UR[3] << CFendl; 
//          Calculates the 2D Euler flux using the exact Riemann Solver
//          -----------------------------------------------------------

//          NB:  Input parameters are the left (UL) and right (UR) state
//          vectors of conserved variables,  with
//          U = (rho, rho u, rho v, rho w, rho E), and the side-area
//          vector, S, where S = d(s1,s2), where d is
//          the area of the cell face and (s1,s2) is the outward
//          pointing cell-face normal vector, pointing from cell L
//          to cell R.
//          The solution is rotated to the S vector frame of reference.
//          The EXACTRS subroutine also needs the xphi and PRESSURE subroutines
//          listed below.
//          The flux is returned in vector US.

    CFint I;
    CFreal SM,ut,vt,wt,uu,vv,ww,sv,pi,qi,ri,RSTAR,CSTAR,EEXP,SEXP,ei;

    CFreal RGAMMA=1.4;
    CFuint conv = 0;
    CFreal eps1 = 1.0E-6;
    CFreal eps2 = 1.0E-6;
    CFreal alpha = 1.0E0;
    CFreal d = std::sqrt((S[2]*S[2]+S[1]*S[1]+S[0]*S[0]));
    CFreal rdd = 1.0E0/(d+1.0E-20);
    CFreal s1 = S[0]*rdd;
    CFreal s2 = S[1]*rdd;
    CFreal s3 = S[2]*rdd;
    CFreal Rrl = 1.0E0/UL[0];
    CFreal Rrr = 1.0E0/UR[0];
    CFreal ql = (s1*UL[1]+s2*UL[2]+s3*UL[3])*Rrl;
    CFreal qr = (s1*UR[1]+s2*UR[2]+s3*UR[3])*Rrr;

//     print*,'##',s1,s2,Rrl,Rrr,ql,qr;

   CFreal pl = pressure(UL);
   CFreal pr = pressure(UR);
   CFreal cl = std::sqrt((RGAMMA*pl*Rrl));
   CFreal cr = std::sqrt((RGAMMA*pr*Rrr));
   CFreal ps = (pl+pr)*0.5;

// CFout << " PL " << pl << " PR " << UL[1] << " CL " << UL[2] << " CR " << UL[3] << CFendl; 

//     print*,'##',pl,pr,cl,cr,ps,ql;

//     !    START OF EXACT RIEMANN SOLVER
   CFreal rconst = UR[0]*std::sqrt(pr);
   CFreal lconst = UL[0]*std::sqrt(pl);
   CFreal MR = rconst*xphi(ps/pr);
   CFreal ML = lconst*xphi(ps/pl);

//     !     Godunov's iteration
  while (conv == 0)
  {
    I=0;
    while ((I <20) && (conv == 0))
    {
      CFreal pn = (pr/MR+pl/ML+ql-qr)/(1.0E0/ML+1.0E0/MR);
      pn = alpha*std::max(eps1,pn)+(1.0E0-alpha)*ps;
      CFreal newmr = rconst*xphi(pn/pr);
      CFreal newml = lconst*xphi(pn/pl);
      CFreal t1 = CFreal(std::abs(newmr-MR));
      CFreal t2 = CFreal(std::abs(newml-ML));
      if ((max(t1,t2) < eps1) || (alpha < eps2)) conv=1;
      ps = pn;
      MR = newmr;
      ML = newml;
      I = I+1;
    }
    alpha = alpha*0.8;
  }
  SM = (pl-pr+ML*ql+MR*qr)/(ML+MR);

//     !print*,'##@!!',SM,pl,pr,ML,ql,MR,qr

//     !    Now find what is happening at the interface
  if (SM < 0.0)
  {
//        !      RIGHT OF THE CONTACT
    ut = UR[1]*Rrr-qr*s1;
    vt = UR[2]*Rrr-qr*s2;
    wt = UR[3]*Rrr-qr*s3;
    if ((ps/pr) > 1.0E0)
    {
//           !      RIGHT TRAVELLING SHOCK
      sv = qr+MR*Rrr;
      if (sv < 0.0E0)
      {
//              !         SUPERSONIC FROM RIGHT TO LEFT
        ri = UR[0];
        qi = qr;
        pi = pr;
      }
      else
      {
//              !         BETWEEN CONTACT AND RIGHT MOVING SHOCK
        ri = MR/(sv-SM);
        qi = SM;
        pi = ps;
      }
    }
    else
    {
//           !      RIGHT TRAVELLING EXPANSION
      EEXP=(qr+cr);
      if (EEXP < 0.0)
      {
//              !         SUPERSONIC FROM RIGHT TO LEFT
        ri = UR[0];
        qi = qr;
        pi = pr;
      }
      else
      {
        CFreal RSTAR=pow((ps*(std::pow(UR[0],RGAMMA))/pr),(1.0E0/RGAMMA));
        CFreal CSTAR=RGAMMA*ps/sqrt(RSTAR);
        CFreal SEXP=(CSTAR+SM);
        if (SEXP > 0.0E0)
        {
//                 !            BETWEEN CONTACT AND START OF EXPANSION
          ri = RSTAR;
          qi = SM;
          pi = ps;
        }
        else
        {
//                 !            IN THE EXPANSION WAVE
          qi=(((RGAMMA-1.0E0)*qr*0.5)-cr)* 2.0E0/(RGAMMA+1.0E0);
          ri=1.0E0+(RGAMMA-1.0E0)*(qi-qr)/(2.0E0*cr);
          ri=UR[0]*ri*(2.0E0/(RGAMMA-1.0E0));
          pi=(pr/(pow(UR[0],RGAMMA)))*(pow(ri,RGAMMA));
        }
      }
    }
  }
  else
  {
//        !      LEFT OF THE CONTACT
    ut = UL[1]*Rrl-ql*s1;
    vt = UL[2]*Rrl-ql*s2;
    wt = UL[3]*Rrl-ql*s3;
    if (ps/pl > 1.0E0)
    {
//           !        LEFT TRAVELLING SHOCK
      sv = ql-ML*Rrl;
      if (sv > 0.0E0)
      {
//              !          SUPERSONIC FROM LEFT TO RIGHT
        qi = ql;
        ri = UL[0];
        pi = pl;
      }
      else
      {
//              !          BETWEEN SHOCK AND CONTACT
        pi = ps;
        qi = SM;
        ri =  ML/(SM-sv);
      }
    }
    else
    {
//           !       LEFT TRAVELLING EXPANSION
      SEXP=ql-cl;
      if (SEXP > 0.0E0)
      {
//              !          Supersonic from left to right
        ri = UL[0];
        pi = pl;
        qi = ql;
      }
      else
      {
        RSTAR=std::pow(ps*( std::pow(UL[0] , RGAMMA) )/pl , 1.0/RGAMMA);
        CSTAR=RGAMMA*ps/std::sqrt(RSTAR);
        EEXP=SM-CSTAR;
//             !       if(ielem .eq. 79) write (*,*) 'SEXP=',SEXP,EEXP
        if (EEXP > 0.0E0)
        {
//                 !             In the expansion wave
          qi=((RGAMMA-1.0E0)*ql*0.5+cl)*2.0E0/(RGAMMA+1.0E0);
          ri=1.0E0-((RGAMMA-1.0E0)/(2.0E0*cl))*(qi-ql);
          ri=UL[0]*std::pow(ri,(2.0E0/(RGAMMA-1.0E0)));
          pi=(pl/(std::pow(UL[0],RGAMMA)))*(pow(ri,RGAMMA));
        }
        else
        {
//                 !             BETWEEN EXPANSION AND CONTACT
          ri = RSTAR;
          qi = SM;
          pi = ps;
//           !    if(ielem .eq. 79) write (*,*) 'ri,qi,pi=',ri,qi,pi
        }
      }
    }
  }


  uu = s1*qi+ut;
  vv = s2*qi+vt;
  ww = s3*qi+wt;

//     !     if(ielem .eq. 79) write (*,*) ri,uu,vv,pi,qi,'&&&&&'
//     !     Find the flux at the interface and scale by the face area, d
  (*US)[0] = ri;
  (*US)[1] = ri*uu;
  (*US)[2] = ri*vv;
  (*US)[3] = ri*ww;
  ei = pi/(RGAMMA-1.0E0)+0.5*ri*(uu*uu+vv*vv+ww*ww);
  (*US)[4] = ei;

}
//////////////////////////////////////////////////////////////////////////////

CFreal ViscousInletBC::xphi(CFreal x)
{
  if (x < 0.999)
  {
   return(kappa1*(1.0-x)/(2.3664319*(1.0- std::pow(x,0.14285714))));
  }
  else
  {
    return(1.2*x+std::sqrt(0.2));
  }
}

//////////////////////////////////////////////////////////////////////////////

CFreal ViscousInletBC::pressure(Framework::State U)
{
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  CFreal pressure;
  if (nbDim == 2)
  {
    pressure = kappa1*(U[3]-(U[1]*U[1]+U[2]*U[2] )/2.0E0/U[0]);
  }
  else
  {
    pressure = kappa1*(U[4]-(U[1]*U[1]+U[2]*U[2]+U[3]*U[3] )/2.0E0/U[0]);
  }
  return pressure;
}

//////////////////////////////////////////////////////////////////////////////

void ViscousInletBC::executeOnTrs()
{
  CFAUTOTRACE;
  CFout << "InletBC applied to " << getCurrentTRS()->getName() << CFendl;
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  RealMatrix T(nbEqs,nbEqs);
  RealMatrix T1(nbEqs,nbEqs);
  RealMatrix Pplus(nbEqs,nbEqs);
  RealMatrix Pminus(nbEqs,nbEqs);
  RealMatrix EigenVal(2,nbEqs);
  RealVector normal;
  State  state;
  State  Dstate;
  State  bstate;
  State  Bstate;
  Node   Dnode;
  CFreal cv = 1.0;

  SafePtr<LinearSystemSolver> lss = getMethodData().getLinearSystemSolver()[0];
  SafePtr<LSSMatrix> jacobMatrix = lss->getMatrix();

  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  Common::SafePtr<GeometricEntityPool<Framework::FaceToCellGEBuilder> >
     geoBuilder = getMethodData().getFaceBuilder();

  // get InnerCells TopologicalRegionSet
  SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells");

  // get Inlet Boundary Faces = CurrentTopologicalRegionSet
  SafePtr<TopologicalRegionSet> faces = getCurrentTRS();

  // get the geodata of the face builder and set the TRSs
  Framework::FaceToCellGEBuilder::GeoData& geoData = geoBuilder->getDataGE();
  geoData.cellsTRS = cells;
  geoData.facesTRS = faces;
  geoData.isBoundary = true;

  //get number of inlet faces
  const CFuint nbFaces = faces->getLocalNbGeoEnts();
//   CFout << "  Faces:" << nbFaces << CFendl;

  //loop over all inlet boundary faces
  for (CFuint iFace = 0; iFace < nbFaces; ++iFace)
  {
    CFLogDebugMax("Face " << iFace << "\n");
    //set index of face
    geoData.idx = iFace;
    //geo builder make face
    m_face = geoBuilder->buildGE();
    //get the (neighbouring) cell
    GeometricEntity* cellLeft  = m_face->getNeighborGeo(LEFT);
    //get nodes of face
    std::vector<Node*>&  nodes  = *m_face->getNodes();
    //get nodes and states of the cell
    std::vector<Node*>&  left_cell_nodes  = *cellLeft->getNodes();
    std::vector<State*>& left_cell_states = *cellLeft->getStates();

    if (left_cell_states[0]->isParUpdatable())
    {
      CFreal detJacobi;
      if (nbDim == 2)
      {
        RealVector hlp = *nodes[1] - *nodes[0];
        detJacobi = sqrt(hlp[0]*hlp[0]+hlp[1]*hlp[1]);
      }
      else
      {
        //     detJacobi = face.computeVolume()*2;
        RealVector hlp1 = *nodes[1] - *nodes[0];
        RealVector hlp2 = *nodes[2] - *nodes[0];
        detJacobi = sqrt((hlp1[1]*hlp2[2] - hlp1[2]*hlp2[1])*(hlp1[1]*hlp2[2] - hlp1[2]*hlp2[1]) + (hlp1[2]*hlp2[0] - hlp1[0]*hlp2[2])*(hlp1[2]*hlp2[0] - hlp1[0]*hlp2[2])+(hlp1[0]*hlp2[1] - hlp1[1]*hlp2[0])*(hlp1[0]*hlp2[1] - hlp1[1]*hlp2[0]));//*2.0/2.0;
      }

      //get number of states of the cell
      const CFuint nbStatesInCellLeft = left_cell_states.size();

      //*****************************************************************
      //*****************************************************************
      //SET Integrator
      //compute shape function in quadrature points
      const std::vector<RealVector>& leftShapeFunctions =  getMethodData().getContourIntegrator()->getSolutionIntegrator(cellLeft)->computeShapeFunctionsAtQuadraturePoints();

      //numbers of quadrature points
      CFuint m_nbKvadrPoint = getMethodData().getContourIntegrator()->getSolutionIntegrator(cellLeft)->getIntegratorPattern()[0];

      //set weights for element quadrature
      const std::vector<RealVector>& leftWeight = getMethodData().getContourIntegrator()->getSolutionIntegrator(cellLeft)->getCoeff();

      //get coordinates of quadrature points
      const std::vector<RealVector>& leftCoord = getMethodData().getContourIntegrator()->getSolutionIntegrator(cellLeft)->getQuadraturePointsCoordinates();

      //compute gradient of shape functions in quadrature points
      std::vector<RealMatrix> leftGradient = cellLeft->computeSolutionShapeFunctionGradientsInMappedCoordinates(leftCoord);

      //*****************************************************************
      //*****************************************************************

      /// @todo this must be improved. Finding in a map
      /// for every cell is efficiently speaking not acceptable.
      //make block acumulator with size of states of the cell
      DGElemTypeData elemData = m_mapElemData.find(nbStatesInCellLeft);
      BlockAccumulator& acc = *elemData.first;
      RealMatrix& elemMat = *elemData.second;
      RealVector& elemVec = *elemData.third;

      //set matrix in blockaccumulator to 0
      elemMat=0.0;
      acc.setValuesM(elemMat);

      // set the IDs on the blockaccumulator (use setRowColIndex() )
      for (CFuint iState = 0; iState < nbStatesInCellLeft; ++iState) {
        const CFuint stateID = left_cell_states[iState]->getLocalID();
        acc.setRowColIndex(iState, stateID);
      }

      //find local index of face in cell (must be improved)
      CFuint m_idxFaceFromLeftCell= 10;
      if (nbDim == 2)
      {
        for (CFuint i=0; i < 3; i++)
          for (CFuint j=0; j < 2; j++)
            if((*left_cell_nodes[(i)%3]==*nodes[(j)%2])&&(*left_cell_nodes[(i+1)%3]==*nodes[(j+1)%2])) m_idxFaceFromLeftCell=i;
      }
      else
      {
        for (CFuint i=0; i < 4; i++)
          for (CFuint j=0; j < 3; j++)
            if(((*left_cell_nodes[(i)%4]==*nodes[(j)%3])&&(*left_cell_nodes[(i+1)%4]==*nodes[(j+1)%3])&&(*left_cell_nodes[(i+2)%4]==*nodes[(j+2)%3]))||((*left_cell_nodes[(i)%4]==*nodes[(j+2)%3])&&(*left_cell_nodes[(i+1)%4]==*nodes[(j+1)%3])&&(*left_cell_nodes[(i+2)%4]==*nodes[(j)%3])))
            {
              m_idxFaceFromLeftCell=i;
            }
      }
      cf_assert(m_idxFaceFromLeftCell!=10);

      //computation of normal (only for straight line connecting boundary nodes of face)
      normal.resize(nbDim);
      if (nbDim == 2)
      {
        normal[0]= (*left_cell_nodes[(m_idxFaceFromLeftCell+1) % 3])[YY] - (*left_cell_nodes[m_idxFaceFromLeftCell])[YY];
        normal[1]= (*left_cell_nodes[m_idxFaceFromLeftCell])[XX] - (*left_cell_nodes[(m_idxFaceFromLeftCell+1) % 3])[XX];
        normal /=sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
      }
      else
      {
        normal = m_face->computeAvgCellNormal();
        RealVector HlpNormal = *left_cell_nodes[(m_idxFaceFromLeftCell+3)%4] - *left_cell_nodes[(m_idxFaceFromLeftCell)%4];
        if ((normal[0]*HlpNormal[0]+normal[1]*HlpNormal[1]+normal[2]*HlpNormal[2]) > 0)
        {
          normal *=-1;
        }
        normal /=sqrt(normal[0]*normal[0] + normal[1]*normal[1] + normal[2]*normal[2]);
      }
      //loop over kvadrature point on the face
      for(CFuint kvadrature_point = 0; kvadrature_point < m_nbKvadrPoint; kvadrature_point++ )
      {
        CFuint leftIndex  = m_idxFaceFromLeftCell*m_nbKvadrPoint + kvadrature_point;
        elemMat=0.0;
        for (CFuint iEq = 0; iEq < nbEqs; ++iEq) //loop over members of state - set to zero
        {
          state[iEq]= 0.0;
          Bstate[iEq]= 0.0;
        }
        //computation of state in point of kvadrature - from previous time step
        for (CFuint iState = 0; iState < nbStatesInCellLeft; ++iState) //loop over states in cell
        {
          RealVector &states = *left_cell_states[iState]->getData();
          for (CFuint iEq = 0; iEq < nbEqs; ++iEq) //loop over members of state
          {
            state[iEq] += leftShapeFunctions[leftIndex][iState]*states[iEq];
          }
        }
        //set coordinates of nodes for gauss kvadrature - to evaluate function from boundary condition
        if (nbDim == 2)
        {
          Dnode = leftCoord[leftIndex][XX]*(*left_cell_nodes[0])
               +leftCoord[leftIndex][YY]*(*left_cell_nodes[1])
               +(1-leftCoord[leftIndex][XX]-leftCoord[leftIndex][YY])*(*left_cell_nodes[2]);
        }
        else
        {
          Dnode = leftCoord[leftIndex][XX]*(*left_cell_nodes[0])
               +leftCoord[leftIndex][YY]*(*left_cell_nodes[1])
               +leftCoord[leftIndex][ZZ]*(*left_cell_nodes[2])
               +(1-leftCoord[leftIndex][XX]   -leftCoord[leftIndex][YY] -leftCoord[leftIndex][ZZ])*(*left_cell_nodes[3]);
        }
        //Dstate = evaluate function from boundary condition in Dnode coordinates
        _vFunction.evaluate(Dnode,Dstate);

        //compute matrixes P+ a P- in point of kvadrature
        if (nbDim == 2)
        {
          compute_EigenValVec2D(state, T, T1, &Pplus, &Pminus, &EigenVal, normal);
        }
        else
        {
          if (compute_EigenValVec3D(state, T, T1, &Pplus, &Pminus, &EigenVal, normal) !=0)
          {
            Node Dnode;
            Dnode = leftCoord[leftIndex][0]*(*(left_cell_nodes[0]))
                 +leftCoord[leftIndex][1]*(*(left_cell_nodes[1]))
                 +leftCoord[leftIndex][2]*(*(left_cell_nodes[2]))
                 +(1-leftCoord[leftIndex][0]-leftCoord[leftIndex][1] -leftCoord[leftIndex][2])*(*(left_cell_nodes[3]));
            CFout << "\n  Negative pressure in " << Dnode << "  STATE  "  << state << "  normal  " << normal << CFendl;
          }
        }
        //compute matrixes K_{sk}
        if (nbDim == 2)
        {
          compute_Kmatrix2D(state,&m_kMatrix);
        }
        else
        {
          compute_Kmatrix3D(state,&m_kMatrix);
        }
        Bstate[0] = state[0];
        for (CFuint iEq = 1; iEq < nbEqs-1; ++iEq) //loop over members of state
        {
          Bstate[iEq] = state[iEq];
          Bstate[nbEqs-1] -= state[iEq]*state[iEq]/(state[0]);
        }
        if (nbDim ==3)
        {
          Bstate[nbEqs-1] = (Bstate[4]/2.0 + state[4])/cv + 1./2.*(Dstate[1]*Dstate[1]+Dstate[2]*Dstate[2]+Dstate[3]*Dstate[3])/Dstate[0];
        }
        else
        {
          Bstate[nbEqs-1] = (Bstate[3]/2.0 + state[3])/cv + 1./2.*(Dstate[1]*Dstate[1]+Dstate[2]*Dstate[2])/Dstate[0];
        }

// if (iFace == 1)  CFout << "\n" << Bstate << "\n";

        //compute inner boundary face term 
        //loop over test function
        for(CFuint row = 0; row < nbStatesInCellLeft; row++ )
        {
          //loop over base function of solution
          for(CFuint col = 0; col < nbStatesInCellLeft; col++ )
          {
            //temp_value = multiplication of test function and base function in kvadrature point
            CFreal temp_value=leftShapeFunctions[leftIndex][row]*leftShapeFunctions[leftIndex][col];
            for(CFuint i = 0; i < nbEqs; i++ )
              for(CFuint j = 0; j < nbEqs; j++ )
              {
                elemMat(row*nbEqs + i, col*nbEqs + j)+=(Pplus(i,j))*temp_value;
              }
          }
        }
        elemMat*=leftWeight[0][kvadrature_point]*detJacobi;
        acc.addValuesM(elemMat);

        elemVec = 0.0;
        //compute stateBC from exact rieman solver (bstate)
        if (nbDim == 2)
        {
          ExactRiemannSolver2D(state, Dstate, normal, &bstate);
        }
        else
        {
          ExactRiemannSolver3D(state, Dstate, normal, &bstate);
        }
        //compute inner boundary face term -> to RHS
        //loop over test function
        for(CFuint row = 0; row < nbStatesInCellLeft; row++ )
        {
          for(CFuint i = 0; i < nbEqs; i++ )//loop over member of state - to test function
          {
            elemVec[row*nbEqs+i]=0; //set RHS[index]
            for(CFuint j = 0; j < nbEqs; j++ )//loop over member of stateBC
            {
              elemVec[row*nbEqs+i]+=(Pminus(i,j))*(bstate)[j];
            }
            elemVec[row*nbEqs+i]*=leftShapeFunctions[leftIndex][row];
          }
        } // end of numerical flux
        for (CFuint iState = 0; iState < nbStatesInCellLeft; ++iState) {
          const CFuint stateID = left_cell_states[iState]->getLocalID();
          for (CFuint iEq = 0; iEq < nbEqs; ++iEq) {
            rhs(stateID, iEq, nbEqs) -= elemVec[iState*nbEqs + iEq]*leftWeight[0][kvadrature_point]*detJacobi;
          }
        }

        elemMat=0.0;
//         if (kvadrature_point==0) elemVec=0.0;
        for(CFuint row = 0; row < nbStatesInCellLeft; row++ )
        {
//         loop over base function of solution
          for(CFuint col = 0; col < nbStatesInCellLeft; col++ )
          {
            for(CFuint k=0;k<nbDim;k++)
            {
              for(CFuint s=0;s<nbDim;s++)
              {
//             temp_value = multiplication of test function, base function in kvadrature point and normal
              CFreal temp_value=leftGradient[leftIndex](col,k)*leftShapeFunctions[leftIndex][row]*normal[s];
              for(CFuint i = 0; i < nbEqs; i++ )
                for(CFuint j = 0; j < nbEqs; j++ )
                {
//                 K_{sk} * \frac{(\partial w}{\partial x_k}|_L * \varphi|_L
                  elemMat(row*nbEqs + i, col*nbEqs + j)+=m_kMatrix[s][k](i,j)*temp_value;
//                  m_theta * K_{sk} * \frac{(\partial \varphi}{\partial x_k}|_L * \w|_L
                  elemMat(col*nbEqs + j, row*nbEqs + i)+=m_theta*m_kMatrix[s][k](i,j)*temp_value;
                }
              }
            }
          }
          //add part from local rhs to global rhs
          const CFuint stateID = left_cell_states[row]->getLocalID();
          for(CFuint k=0;k<nbDim;k++)
          {
            for(CFuint s=0;s<nbDim;s++)
            {
              for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
              {
                for (CFuint j = 0; j < nbEqs; ++j)
                {
                  rhs(stateID, iEq, nbEqs) -= m_theta*m_kMatrix[s][k](j,iEq)*leftGradient[leftIndex](row,k)*Bstate[j]*normal[s]*leftWeight[0][kvadrature_point]*detJacobi;
                }
              }
            }
          }
          for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
          {
            rhs(stateID, iEq, nbEqs) += m_sigma*leftShapeFunctions[leftIndex][row]*Bstate[iEq]*leftWeight[0][kvadrature_point]*detJacobi/detJacobi/m_Re;
//             elemVec[row + iEq*nbStatesInCellLeft] += m_sigma*leftShapeFunctions[leftIndex][row]*Bstate[iEq]*leftWeight[0][kvadrature_point]*detJacobi/detJacobi/m_Re;
          }
          for(CFuint col = 0; col < nbStatesInCellLeft; col++ )
          {
            for (CFuint iEq = 0; iEq < nbEqs; ++iEq)
            {
              elemMat(row*nbEqs + iEq, col*nbEqs + iEq)-=m_sigma*leftShapeFunctions[leftIndex][row]*leftShapeFunctions[leftIndex][col]/detJacobi/m_Re;
            }
          }
        }
        elemMat*=-leftWeight[0][kvadrature_point]*detJacobi;
        acc.addValuesM(elemMat);
      }
//       if (iFace == 1)
//          CFout << "\n" << elemVec << "\n";
      // add the values in the jacobian matrix
      jacobMatrix->addValues(acc);
    }
    // release the face
    geoBuilder->releaseGE();
  }
  CFout << " ... OK\n" << CFendl;

//       jacobMatrix->finalAssembly();
//       jacobMatrix->printToFile("inside");

}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> > ViscousInletBC::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace DiscontGalerkin
}  // namespace COOLFluiD

