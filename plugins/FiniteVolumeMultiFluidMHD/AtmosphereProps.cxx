#include "Common/PE.hh"
#include "Common/BadValueException.hh"

#include "MathTools/MathConsts.hh"

#include "Environment/FileHandlerOutput.hh"
#include "Environment/DirPaths.hh"
#include "Environment/SingleBehaviorFactory.hh"

#include "Framework/PathAppender.hh"
#include "Framework/DataProcessing.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalChemicalLibrary.hh"
#include "Framework/PhysicalConsts.hh"

#include "FiniteVolume/CellCenterFVM.hh"
#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/AtmosphereProps.hh"


/////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Numerics::FiniteVolume;
//using namespace COOLFluiD::Physics::Maxwell;
//using namespace COOLFluiD::Physics::MultiFluidMHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<AtmosphereProps,
    DataProcessingData, FiniteVolumeMultiFluidMHDModule>
AtmospherePropsProvider("AtmosphereProps");

//////////////////////////////////////////////////////////////////////////////

void AtmosphereProps::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

AtmosphereProps::AtmosphereProps(const std::string& name) :
  DataProcessingCom(name),
  socket_states("states"),
  socket_gstates("gstates"),
  socket_nstates("nstates"),
  socket_volumes("volumes"),
  socket_nodes("nodes"),
  socket_isOutward("isOutward"),
  socket_normals("normals"),
  socket_gradPi("gradPi"),
  socket_gradPn("gradPn"),
  socket_nuIon("nuIon"),
  socket_nuRec("nuRec"),
  socket_nu_in("nu_in"),
  socket_jxB_x("jxB_x"),
  socket_jxB_y("jxB_y"),
  socket_jxB_z("jxB_z"),
  socket_jxB("jxB"),
  socket_Jx("Jx"),
  socket_Jy("Jy"),
  socket_Jz("Jz"),
  socket_Jtot("Jtot"),
  m_gradPi(),
  m_gradPn(),
  m_curlB(),
  m_normal()

{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

AtmosphereProps::~AtmosphereProps()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
AtmosphereProps::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  
  result.push_back(&socket_gradPi);
  result.push_back(&socket_gradPn);
  result.push_back(&socket_nuIon);
  result.push_back(&socket_nuRec);
  result.push_back(&socket_nu_in);
  result.push_back(&socket_Jx);
  result.push_back(&socket_Jy);
  result.push_back(&socket_Jz);
  result.push_back(&socket_Jtot);
  result.push_back(&socket_jxB_x);
  result.push_back(&socket_jxB_y);
  result.push_back(&socket_jxB_z);
  result.push_back(&socket_jxB);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
AtmosphereProps::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  
  result.push_back(&socket_states);
  result.push_back(&socket_gstates);
  result.push_back(&socket_nstates);
  result.push_back(&socket_volumes);
  result.push_back(&socket_nodes);
  result.push_back(&socket_isOutward);
  result.push_back(&socket_normals);  
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void AtmosphereProps::setup()
{
  CFAUTOTRACE;
  
  // Get number of cells
  Common::SafePtr<Common::ConnectivityTable<CFuint> > cells =
    MeshDataStack::getActive()->getConnectivity("cellStates_InnerCells");

  const CFuint nbCells = cells->nbRows();
  const CFuint DIM = 2; //AAL: to generalize in the future for 3D

  DataHandle<CFreal> gradPi = socket_gradPi.getDataHandle();
  gradPi.resize(nbCells);
  gradPi = 0.0;

  DataHandle<CFreal> gradPn = socket_gradPn.getDataHandle();
  gradPn.resize(nbCells);
  gradPn = 0.0;
  
  DataHandle<CFreal> nuIon = socket_nuIon.getDataHandle();
  nuIon.resize(nbCells);
  nuIon = 0.0;
  
  DataHandle<CFreal> nuRec = socket_nuRec.getDataHandle();
  nuRec.resize(nbCells);
  nuRec = 0.0;
  
  DataHandle<CFreal> nu_in = socket_nu_in.getDataHandle();
  nu_in.resize(nbCells);
  nu_in = 0.0;

  DataHandle<CFreal> Jx = socket_Jx.getDataHandle();
  Jx.resize(nbCells);
  Jx = 0.0;

  DataHandle<CFreal> Jy = socket_Jy.getDataHandle();
  Jy.resize(nbCells);
  Jy = 0.0;

  DataHandle<CFreal> Jz = socket_Jz.getDataHandle();
  Jz.resize(nbCells);
  Jz = 0.0;

  DataHandle<CFreal> Jtot = socket_Jtot.getDataHandle();
  Jtot.resize(nbCells);
  Jtot = 0.0;

  DataHandle<CFreal> jxB_x = socket_jxB_x.getDataHandle();
  jxB_x.resize(nbCells);
  jxB_x = 0.0;

  DataHandle<CFreal> jxB_y = socket_jxB_y.getDataHandle();
  jxB_y.resize(nbCells);
  jxB_y = 0.0;

  DataHandle<CFreal> jxB_z = socket_jxB_z.getDataHandle();
  jxB_z.resize(nbCells);
  jxB_z = 0.0;

  DataHandle<CFreal> jxB = socket_jxB.getDataHandle();
  jxB.resize(nbCells);
  jxB = 0.0;

  m_geoBuilder.setup();

// AAL: To add in the future a chemical library to compute this
//  m_library = PhysicalModelStack::getActive()->getImplementor()->
//    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
//  cf_assert(m_library.isNotNull());
  
  m_gradPi.resize(DIM, 0.);
  m_gradPn.resize(DIM, 0.);
  m_curlB.resize(3, 0.);
  m_normal.resize(DIM, 0.);
}

//////////////////////////////////////////////////////////////////////////////

void AtmosphereProps::execute()
{
  CFout <<"AtmosphereProps::computing properties \n";
  CFAUTOTRACE;
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  DataHandle<RealVector> nstates = socket_nstates.getDataHandle();
  DataHandle<Framework::State*> gstates = socket_gstates.getDataHandle();
  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();
  DataHandle<CFint> isOutward = socket_isOutward.getDataHandle();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  
  DataHandle<CFreal> gradPi = socket_gradPi.getDataHandle();
  DataHandle<CFreal> gradPn = socket_gradPn.getDataHandle();
  DataHandle<CFreal> nuIon = socket_nuIon.getDataHandle();
  DataHandle<CFreal> nuRec = socket_nuRec.getDataHandle();
  DataHandle<CFreal> nu_in = socket_nu_in.getDataHandle();

  DataHandle<CFreal> Jx = socket_Jx.getDataHandle();
  DataHandle<CFreal> Jy = socket_Jy.getDataHandle();
  DataHandle<CFreal> Jz = socket_Jz.getDataHandle();
  DataHandle<CFreal> Jtot = socket_Jtot.getDataHandle();
  DataHandle<CFreal> jxB_x = socket_jxB_x.getDataHandle();
  DataHandle<CFreal> jxB_y = socket_jxB_y.getDataHandle();
  DataHandle<CFreal> jxB_z = socket_jxB_z.getDataHandle();
  DataHandle<CFreal> jxB = socket_jxB.getDataHandle();
  
  const CFuint DIM = 2; //AAL: to generalize to 3D in the future
  
  m_geoBuilder.getGeoBuilder()->setDataSockets(socket_states, socket_gstates, socket_nodes);
  
  std::vector<SafePtr<TopologicalRegionSet> > trsList = MeshDataStack::getActive()->getTrsList();

  Common::SafePtr<TopologicalRegionSet> cells = MeshDataStack::getActive()->getTrs("InnerCells"); 
  CellTrsGeoBuilder::GeoData& geoData = m_geoBuilder.getDataGE();
  geoData.trs = cells;

  const CFuint nbCells = cells->getLocalNbGeoEnts();

  for (CFuint iCell = 0; iCell < nbCells; ++iCell) {
    
    geoData.idx = iCell;
    GeometricEntity* currCell = m_geoBuilder.buildGE();
    const CFuint elemID = currCell->getID();
    m_gradPi = 0.;
    m_gradPn = 0.;
    m_curlB  = 0.;
    
    State *currState = currCell->getState(0);

    CFreal nu_Ion, nu_Rec, ionsIonizRate, neutralsRecombRate;
    CFreal nu_inLoc, nu_niLoc, nu_enLoc, nu_eiLoc;

    computeChemFreqs(currState, nu_Ion, nu_Rec,
                ionsIonizRate, neutralsRecombRate);

    computeCollFreqs(currState, nu_inLoc, nu_niLoc, nu_enLoc, nu_eiLoc);


    const vector<GeometricEntity*>& faces = *currCell->getNeighborGeos();
    const CFuint nbFaces = faces.size();
    const CFreal ovVolume = 1./volumes[elemID];

    for (CFuint i = 0; i < nbFaces; ++i) {
      const GeometricEntity *const face = currCell->getNeighborGeo(i);
      const CFuint faceID = face->getID();
      const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();
      const CFuint nbFaceNodes = face->nbNodes();
      const CFreal ovNbFaceNodes = 1./(CFreal)nbFaceNodes;
      
      // store the outward face normal (scaled with the corresponding face area)
      // fill in only the components that are available (==dim)
      const CFreal factor = (static_cast<CFuint>(isOutward[faceID]) != elemID) ? -1. : 1.;
      for (CFuint d = 0; d < DIM; ++d) {
        m_normal[d] = normals[startID+d]*factor;
      }
      const CFreal nx = m_normal[0];
      const CFreal ny = m_normal[1];
      const CFreal nz = 0; //AAL: Default 2D

      for (CFuint n = 0; n < nbFaceNodes; ++n) {
        const CFuint nodeID = face->getNode(n)->getLocalID();
        const RealVector& nodalState = nstates[nodeID];
        // consider all 3D components here even in 2D
        for (CFuint d = 0; d < DIM; ++d) {
          const CFreal nd = m_normal[d];
          const CFuint rhoiID = 8;
          const CFuint rhonID = 9;
          const CFuint TiID   = 14;
          const CFuint TnID   = 15;

          const CFreal rhoiNode = nodalState[rhoiID];
          const CFreal rhonNode = nodalState[rhonID];
          const CFreal TiNode = nodalState[TiID];
          const CFreal TnNode = nodalState[TnID];

          const CFreal PiNodal = getIonPressure(rhoiNode, TiNode);
          const CFreal PnNodal = getNeutralPressure(rhonNode, TnNode);

          m_gradPi[d] += nd*PiNodal*ovNbFaceNodes;
          m_gradPn[d] += nd*PnNodal*ovNbFaceNodes;  
        } //Loop over the dimension
        const CFuint BxID = 0;
        const CFuint ByID = 1;
        const CFuint BzID = 2;

        const CFreal BxNode = nodalState[BxID];
        const CFreal ByNode = nodalState[ByID];
        const CFreal BzNode = nodalState[BzID];

        m_curlB[0] += (-nz*ByNode + ny*BzNode)*ovNbFaceNodes;
        m_curlB[1] += (-nx*BzNode + nz*BxNode)*ovNbFaceNodes;
        m_curlB[2] += (-ny*BxNode + nx*ByNode)*ovNbFaceNodes;

      } //Loop over the nodes
    } //Loop over the faces

    m_gradPi *= ovVolume;
    m_gradPn *= ovVolume;
    m_curlB  *= ovVolume;
    
    gradPi[iCell] = m_gradPi[1];
    gradPn[iCell] = m_gradPn[1];

    const CFreal ovmu0 = 1./1.2566370614e-6; //mu0 is the magnetic permeability in units [H/m], source Wikipedia
    const CFreal Bx = (*currState)[0];
    const CFreal By = (*currState)[1];
    const CFreal Bz = (*currState)[2];

    Jx[iCell]   = m_curlB[0]*ovmu0;
    Jy[iCell]   = m_curlB[1]*ovmu0;
    Jz[iCell]   = m_curlB[2]*ovmu0;
    Jtot[iCell] = std::sqrt(Jx[iCell]*Jx[iCell] + Jy[iCell]*Jy[iCell]
                + Jz[iCell]*Jz[iCell]);

    jxB_x[iCell] = Jy[iCell]*Bz - Jz[iCell]*By;
    jxB_y[iCell] = Jz[iCell]*Bx - Jx[iCell]*Bz;
    jxB_z[iCell] = Jx[iCell]*By - Jy[iCell]*Bx;
    jxB[iCell]   = std::sqrt(jxB_x[iCell]*jxB_x[iCell] + jxB_y[iCell]*jxB_y[iCell]
                + jxB_z[iCell]*jxB_z[iCell]);

    nuIon[iCell] = nu_Ion;
    nuRec[iCell] = nu_Rec;
    nu_in[iCell] = nu_inLoc;
    
    // Set the boundary at 0
//    for (CFuint iFace = 0; iFace < nbFaces; ++iFace) {
//      const GeometricEntity *const face = currCell->getNeighborGeo(iFace);
//      State *const neighborState = (currCell->getState(0) == face->getState(0)) ? face->getState(1) : face->getState(0);
//      if (neighborState->isGhost()){
//        //gradPi[iCell] = 0;
//        //gradPn[iCell] = 0;
//      }
//    }
    m_geoBuilder.releaseGE();
  }
}
//////////////////////////////////////////////////////////////////////////////
void AtmosphereProps::getDensTemp(const Framework::State* currState, CFreal& rhoi, CFreal& rhon, CFreal& Ti, CFreal& Tn)
{
    const CFint endEM = 8;
    const CFint DIM = 2;
    const CFint nbSpecies = 2;
    rhoi = (*currState)[endEM];
    rhon = (*currState)[endEM + 1];
    Ti   = (*currState)[endEM + nbSpecies + DIM*nbSpecies];
    Tn   = (*currState)[endEM + nbSpecies + DIM*nbSpecies + 1];
}

//////////////////////////////////////////////////////////////////////////////

void AtmosphereProps::computeChemFreqs(const Framework::State* currState,
        CFreal& nu_Ion, CFreal& nu_Rec, CFreal& ionsIonizRate, CFreal& neutralsRecombRate)
{
  //cout <<"AtmosphereProps::computeChemFreqs \n";
  CFreal rhoi, rhon, Ti, Tn;

  getDensTemp(currState, rhoi, rhon, Ti, Tn);

  //Molecular Masses
  //const CFreal me = 9.1094e-31;              // Electron's mass [kg] source:Standart Handbook for Electrical Engineerings
  const CFreal mi = 1.6726e-27;              // Proton's mass [kg] source:Standart Handbook for Electrical Engineerings
  const CFreal mn = 1.6726e-27;              // Neutral's mass [kg] source:Standart Handbook for Electrical Engineerings

  //electron properties
  const CFreal ne = rhoi/mi;		   	// Electrically neutral, i.e. ne = ni
  const CFreal Te = Ti/11604.50520;		// electrons temperature in eV. Thermal equilibrium is assumed Ti = Te

  //Neutrals and ions properties
  const CFreal nn = rhon/mn;			// neutral particle per unit volume
  const CFreal ni = ne;			        // ion particle per unit volume
  // IONIZATION
  //constants Related to ionization frequency from [Leake]

  nu_Ion = 0.;
  //bool correctedChemistry = true; //TODO: To be added as option

  //if(correctedChemistry) {
    const CFreal factor = 0.6;
    const CFreal beta = 158000*factor/Ti;
    const CFreal I = 2.34e-14*std::pow(beta,-0.5)*std::exp(-beta);

    nu_Ion = ne*I;

//  }
//  else {
//    const CFreal A = 2.91e-14;
//    const CFreal X = 0.232;
//    const CFreal psiIonOvTe = 13.6/Te;
//    const CFreal K = 0.39;

//    nu_Ion = ne*A/(X + psiIonOvTe)*std::pow(psiIonOvTe, K)*std::exp(-psiIonOvTe);  // Ionization freq.
//  }
  const CFreal GammaIon_n = -nn*nu_Ion;

  // RECOMBINATION
  nu_Rec = 0.;
  //constant related to recombination
//  if(correctedChemistry) {
//    const CFreal factor = 0.6;
//    const CFreal beta = 158000*factor/Te;
    const CFreal R = 5.20e-20*std::pow(beta,0.5)*(0.4288 + 0.5*std::log(beta) + 0.4698*std::pow(beta,-1/3));

    nu_Rec = ne*R;
//  }
//  else {
//    const CFreal B = 2.6e-19;
//    nu_Rec = ne/std::sqrt(Te)*B;
//  }
 const CFreal GammaRec_i = -ni*nu_Rec;

  ///TOTAL (particles/m3)
  ionsIonizRate     = -GammaIon_n;
  neutralsRecombRate = -GammaRec_i;
}

//////////////////////////////////////////////////////////////////////////////

void AtmosphereProps::computeCollFreqs(const Framework::State* currState,
        CFreal& nu_in, CFreal& nu_ni, CFreal& nu_en, CFreal& nu_ei)
{
  //cout <<"AtmosphereProps::computeCollFreqs \n";
  // particles density, temperature
  CFreal rhoi, rhon, Ti, Tn;

  getDensTemp(currState, rhoi, rhon, Ti, Tn);

  //data
  // particle mass
  const CFreal mi = 1.6726e-27;  // Proton's mass [kg] source:Standart Handbook for Electrical Engineerings
  const CFreal mn = 1.6726e-27;  // Neutral's mass [kg] source:Standart Handbook for Electrical Engineerings
  const CFreal me = 9.1093829140e-31;
  const CFreal kB = Framework::PhysicalConsts::Boltzmann(); // Boltzmann constant
  const CFreal e_charge = 1.60217656535e-19; //Electron's charge
  const CFreal epsilon  = 8.8542e-12;

  // particles per volume
  const CFreal ni = rhoi/mi;
  const CFreal nn = rhon/mn;
  const CFreal ne = ni;

  // electron-neutral collision frequency
  const CFreal Sigma_en   = 1e-19;
  const CFreal pi         = MathTools::MathConsts::CFrealPi(); //Pi number
  const CFreal m_en       = me*mn/(me + mn);
  const CFreal T_en       = (Ti + Tn)/2;
  nu_en      = nn*Sigma_en*std::sqrt(8.*kB*T_en/(pi*m_en));

  // electron-ion collision frequency
  const CFreal T_ei       = Ti; //assumed thermal equilibrium
  const CFreal r_debye    = e_charge*e_charge/(4*pi*epsilon*kB*Ti);
  const CFreal Sigma_ei   = pi*r_debye*r_debye;
  const CFreal m_ei       = me*mi/(me + mi);
  nu_ei      = ni*Sigma_ei*std::sqrt(8.*kB*T_ei/(pi*m_ei));

  //parameters
  const CFreal m_in = mi*mn/(mi + mn);
  const CFreal Sigma_in = 1.41e-19;		//collisional cross-section m2 [Leake]
  const CFreal T_in = (Ti + Tn)/2;

  //ion-neutral collision Frequency
  nu_in = nn*Sigma_in*std::sqrt(8.*kB*T_in/(pi*m_in));
  nu_ni = nu_in/nn*ni;

}

//////////////////////////////////////////////////////////////////////////////

void AtmosphereProps::unsetup()
{
  CFAUTOTRACE;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

