
#include "Environment/DirPaths.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/SystemMITReM.hh"
#include "Muffin/SystemMITReMPLaS.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< SystemMITReMPLaS,MuffinData,MuffinModule > cSystemMITReMPLaSProvider("MITReMPLaS");

//////////////////////////////////////////////////////////////////////////////

SystemMITReMPLaS::SystemMITReMPLaS(const std::string& name) :
    SystemPLaS(name),
    s_gasonsurf("GasOnSurface")  // socket sinks
{
  CFAUTOTRACE;

  // configuration options for this command (calls defineConfigOptions)
  attachTag("Muffin::SystemMITReMPLaS");
  addConfigOptionsTo(this);

  // gas evolution
  m_ge_fmax = 0.9;
  m_bubbles_type_str.clear();
  m_bubbles_mu  = 0.;
  m_bubbles_sig = 0.;
  m_bubbles_min = 0.;
  m_bubbles_t   = 298.15;
  setParameter("GEvolutionSurfaceFractionMax",&m_ge_fmax);
  setParameter("BubblesDType",&m_bubbles_type_str);
  setParameter("BubblesDMean",&m_bubbles_mu);
  setParameter("BubblesDStdDev",&m_bubbles_sig);
  setParameter("BubblesDMinDiam",&m_bubbles_min);
  setParameter("BubblesGTemperature",&m_bubbles_t);
}

//////////////////////////////////////////////////////////////////////////////

void SystemMITReMPLaS::defineConfigOptions(Config::OptionList& options)
{
  // gas evolution
  options.addConfigOption< double >("GEvolutionSurfaceFractionMax","Gas evolution surface gas fraction maximum (from 0. to 1., default 0.9)");
  options.addConfigOption< std::string >("BubblesDType","Bubbles distribution type (\"None\" (default), \"Constant\", \"Normal\" or \"Log-Normal\")");
  options.addConfigOption< double >("BubblesDMean","Bubbles distribution mean (default 0.)");
  options.addConfigOption< double >("BubblesDStdDev","Bubbles distribution standard deviation (default 0.)");
  options.addConfigOption< double >("BubblesDMinDiam","Bubbles distribution minimum diameter (default 0.)");
  options.addConfigOption< double >("BubblesGTemperature","Generated bubbles temperature (default 298.15)");
}

//////////////////////////////////////////////////////////////////////////////

void SystemMITReMPLaS::setup()
{
  CFAUTOTRACE;
  SystemPLaS::setup();


  // inportant DataHandle's
  DataHandle< std::vector< GasOnSurface > > h_gasonsurf = s_gasonsurf.getDataHandle();
  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  const CFuint nbDim = h_nodes[0]->size();

  std::vector< SafePtr< TopologicalRegionSet > > btrs = MeshDataStack::getActive()->getFilteredTrsList("boundary");


  log("set gas evolution properties...");
  m_ge_fmax = std::max(0.,std::min(1.,m_ge_fmax));
  m_bubbles_type = (m_bubbles_type_str=="Constant"?   1:
                   (m_bubbles_type_str=="Normal"?     2:
                   (m_bubbles_type_str=="Log-Normal"? 3:0 )));
  if (m_bubbles_type_str=="Log-Normal" && m_bubbles_mu<=MathTools::MathConsts::CFrealEps())
    err("for BubblesDType=Log-Normal BubblesDMean must be > 0.");
  if (m_bubbles_type_str=="Log-Normal" && m_bubbles_sig<=MathTools::MathConsts::CFrealEps())
    err("for BubblesDType=Log-Normal BubblesDStdDev must be > 0.");
  ver("set gas evolution properties.");


  log("set initial bubbles distribution...");
  m_vbubbles.resize(h_gasonsurf.size());
  for (CFuint i=0; i<h_gasonsurf.size(); ++i) {
    if (h_gasonsurf[i].size()) {
      const ConnectivityTable< CFuint >& faces = *btrs[i]->getGeo2NodesConn();

      const Bubble bempty = { 0.,0.,0., 0.,0.,0., m_bubbles_mu, m_bubbles_t };
      m_vbubbles[i].assign(faces.nbRows(),bempty);

      for (CFuint f=0; f<faces.nbRows(); ++f)
        if (h_gasonsurf[i][f].isLocal)
          m_vbubbles[i][f] = generateNewBubble(
            h_nodes[faces(f,0)],
            h_nodes[faces(f,1)],
            (nbDim>2? h_nodes[faces(f,2)] : CFNULL) );
    }
  }
  ver("set initial bubbles distribution.");
}

//////////////////////////////////////////////////////////////////////////////

void SystemMITReMPLaS::execute()
{
  CFAUTOTRACE;

  // get DataHandle's
  DataHandle< std::vector< GasOnSurface > > h_gasonsurf = s_gasonsurf.getDataHandle();
  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  const CFuint nbDim = h_nodes[0]->size();

  std::vector< SafePtr< TopologicalRegionSet > > btrs = MeshDataStack::getActive()->getFilteredTrsList("boundary");


  // check there is a timestep set
  SafePtr< SubSystemStatus > status = getMethodData().m_status;
  const CFreal dt = status->getDTDim();
  if (dt<MathTools::MathConsts::CFrealEps()) {
    err("dt is 0., Loop is likely not correct");
  }


  log("create bubbles from accumulated gas...");
  std::vector< double > vextEntPos;
  std::vector< double > vextEntVel;
  std::vector< double > vextEntDiam;
  std::vector< double > vextEntTemp;
  for (CFuint b=0; b<btrs.size(); ++b) {
    if (!btrs[b]->hasTag("Muffin::GasProduction"))
      continue;

    const ConnectivityTable< CFuint >& faces = *btrs[b]->getGeo2NodesConn();
    double dVdt = 0.;  // gas  production rate
    double Vsrf = 0.;  //  ... on surface
    double Vrls = 0.;  //  ... released into the flow
    for (CFuint f=0; f<faces.nbRows(); ++f) {
      GasOnSurface& gasonsurf = h_gasonsurf[b][f];
      if (!gasonsurf.isLocal)
        continue;

      // time integration
      gasonsurf.V += gasonsurf.dVdt * dt;

      // generate bubbles while accumulated gas volume bigger than
      // to-be-bubble volume
      double Vdetach = M_PI/6.*pow(m_vbubbles[b][f].d,3.);
      while (gasonsurf.V>Vdetach) {

        // detach old bubble
        Bubble& bubble = m_vbubbles[b][f];
        vextEntPos.push_back(bubble.x);
        vextEntPos.push_back(bubble.y);
        vextEntPos.push_back(bubble.z);
        vextEntVel.push_back(bubble.u);
        vextEntVel.push_back(bubble.v);
        vextEntVel.push_back(bubble.w);
        vextEntDiam.push_back(bubble.d);
        vextEntTemp.push_back(bubble.t);
        gasonsurf.V -= Vdetach;
        Vrls        += Vdetach;

        // create new bubble
        m_vbubbles[b][f] = generateNewBubble(
          h_nodes[faces(f,0)],
          h_nodes[faces(f,1)],
          (nbDim>2? h_nodes[faces(f,2)] : CFNULL) );
        Vdetach = M_PI/6.*pow(m_vbubbles[b][f].d,3.);

      }

      // update surface gas fraction
      const double S = M_PI/4. * pow( 6.*gasonsurf.V/M_PI , 2./3. );
      gasonsurf.F = std::min(S/gasonsurf.S,m_ge_fmax);

      Vsrf += gasonsurf.V;     // update gas on surface
      dVdt += gasonsurf.dVdt;  // ... gas production rate
    }

    log( "TRS " + btrs[b]->getName() + " gas: "
      + " p. rate [m3 s-1] / v. surface [m3] / v. released [m3]: "
      + StringOps::to_str(dVdt) + " / "
      + StringOps::to_str(Vsrf) + " / "
      + StringOps::to_str(Vrls) );

  }
  log("number: " + StringOps::to_str(vextEntDiam.size()));
  ver("create bubbles from accumulated gas.");


  log("write surface gas tracking/bubbles information...");
  writeAttachedBubbles(
    m_plas->getInputParameters().writeTecplotFilename,
    status->getNbIter(),
    status->getCurrentTimeDim() );
  ver("write surface gas tracking/bubbles information.");


  log("process...");
  const int numExtEnt = (int) vextEntDiam.size();
  m_plas->setIterationProperties(
    status->getNbIter(),
    status->getCurrentTimeDim(),
    status->getDTDim(),
    1, numExtEnt,
    (!numExtEnt? CFNULL:&vextEntPos[0]),
    (!numExtEnt? CFNULL:&vextEntVel[0]),
    (!numExtEnt? CFNULL:&vextEntTemp[0]),
    (!numExtEnt? CFNULL:&vextEntDiam[0]) );
  m_plas->unblock();  // m_plas is mine, all mine!
  m_plas->process();  // do it, yeah
  m_plas->block();    // block intruders
  ver("process.");


  log("update void fraction...");
  DataHandle< CFreal > h_mn_voidfraction = s_mn_voidfraction.getDataHandle();
  const PLAS_PHASE_DATA* m_plas_phase = m_plas->getPhaseData();
  for (CFuint n=0; n<h_mn_voidfraction.size(); ++n) {
    const double vf = m_plas_phase[n].volFrac;
    h_mn_voidfraction[n] = (CFreal) std::min(1.,std::max(0.,vf));
  }
  ver("update void fraction.");
}

//////////////////////////////////////////////////////////////////////////////

double SystemMITReMPLaS::generateRGaussian(double m, double s) const
{
  static bool use_last = false;
  static double yy2;
  double yy1 = yy2;

  if (!use_last) {
    double xx1, xx2, w;
    do {
      xx1 = generateRDouble(-1.,1.);
      xx2 = generateRDouble(-1.,1.);
      w = xx1*xx1 + xx2*xx2;
    }
    while (w>=1.);

    w = sqrt( (-2.*::log(w)) / w );
    yy1 = xx1*w;
    yy2 = xx2*w;
  }

  use_last = !use_last;
  return m + yy1*s;
}

//////////////////////////////////////////////////////////////////////////////

Bubble SystemMITReMPLaS::generateNewBubble(Node* p1, Node* p2, Node* p3)
{
  Bubble b = {
    0.,0.,0.,       // original position
    0.,0.,0.,       // initial velocity
    m_bubbles_mu,   // diameter
    m_bubbles_t };  // temperature

  // generate new position
  const CFuint nbDim = p1->size();
  if (nbDim>2) {
    double r1 =  0.;
    double r2 = 10.;
    while (r1+r2>1.) {
      r1 = generateRDouble(0.,1.);
      r2 = generateRDouble(0.,1.);
    }
    b.x = (*p1)[0] + ((*p2)[0]-(*p1)[0])*r1 + ((*p3)[0]-(*p1)[0])*r2;
    b.y = (*p1)[1] + ((*p2)[1]-(*p1)[1])*r1 + ((*p3)[1]-(*p1)[1])*r2;
    b.z = (*p1)[2] + ((*p2)[2]-(*p1)[2])*r1 + ((*p3)[2]-(*p1)[2])*r2;
  }
  else {
    const double r = generateRDouble(0.,1.);
    b.x = (*p1)[0] + ((*p2)[0]-(*p1)[0])*r;
    b.y = (*p1)[1] + ((*p2)[1]-(*p1)[1])*r;
    b.z = 0.;
  }

  // generate new diameter
  if (m_bubbles_type==2) {  // ... normal spectrum

    do {
      b.d = generateRGaussian(m_bubbles_mu,m_bubbles_sig);
    } while (b.d<m_bubbles_min);

  }
  else if (m_bubbles_type==3) {  // ... log-normal spectrum
    const double sigLn = sqrt(::log(
      (m_bubbles_sig*m_bubbles_sig/exp(2.*::log(m_bubbles_mu))) + 1. ));
    const double muLn = ::log(m_bubbles_mu) - 0.5*sigLn*sigLn;
    do {
      b.d = ::exp(generateRGaussian(muLn,sigLn));
    } while (b.d<m_bubbles_min);

  }

  return b;
}

//////////////////////////////////////////////////////////////////////////////

void SystemMITReMPLaS::writeAttachedBubbles(const std::string& fname, const unsigned i, const double t)
{
  DataHandle< Node*,GLOBAL > h_nodes = s_nodes.getDataHandle();
  DataHandle< std::vector< GasOnSurface > > h_gasonsurf = s_gasonsurf.getDataHandle();
  const CFuint nbDim = h_nodes[0]->size();

  // count number of sticking bubbles
  unsigned numEnt = 0;
  for (unsigned b=0; b<m_vbubbles.size(); ++b)
    numEnt += m_vbubbles[b].size();
  GlobalReduceOperation< GRO_SUM >(&numEnt,&numEnt);

  // write attached bubbles, each rank in turn
  const unsigned numProc = PE::GetPE().GetProcessorCount();
  const unsigned irank   = PE::GetPE().GetRank();
  for (unsigned iproc=0; iproc<numProc; ++iproc) {
    PE::GetPE().setBarrier();
    if (iproc!=irank)
      continue;

    // first rank writes zone header, then all ranks write their bubbles
    std::ofstream o(fname.c_str(),std::ios_base::app);
    o.precision(16);
    if (!irank)
      o << "ZONE T=\"Entities (attached) iter=" << i << " time=" << t << "\""
        << " I=" << numEnt
        << " SOLUTIONTIME=" << t
        << " DATAPACKING=POINT" << std::endl;
    for (unsigned b=0; b<m_vbubbles.size(); ++b) {
      for (unsigned f=0; f<m_vbubbles[b].size(); ++f) {
        const Bubble& B = m_vbubbles[b][f];
        const double  d = pow(h_gasonsurf[b][f].V*6./M_PI,1./3.);
        o << B.x << ' ' << B.y << ' ' << (nbDim>2? StringOps::to_str(B.z):"") << ' '
          << B.u << ' ' << B.v << ' ' << (nbDim>2? StringOps::to_str(B.w):"") << ' '
          << d   << ' ' << B.t << std::endl;
      }
    }
    o.close();
  }
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

