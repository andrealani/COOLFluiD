#include "Framework/MethodStrategyProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/DataHandle.hh"
#include "Framework/GeometricEntity.hh"
#include "MHD/MHDTerm.hh"
#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"
#include "FiniteVolumeMHD/MHD2DConsACAHWSourceTerm.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

////////////////////////////////////////////////////////////////////////////// 

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

////////////////////////////////////////////////////////////////////////////// 

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<MHD2DConsACAHWSourceTerm,
		       CellCenterFVMData,
		       ComputeSourceTerm<CellCenterFVMData>,
		       FiniteVolumeMHDModule>
mHD2DConsACAHWSTFVMCCProvider("MHD2DConsACAHWST");

//////////////////////////////////////////////////////////////////////////////

MHD2DConsACAHWSourceTerm::MHD2DConsACAHWSourceTerm(const std::string& name) :
  ComputeSourceTermFVMCC(name),
  _gradP(),
  _gradBx(),
  _gradBy(),
  //_gradBz(),
  _gradRho(),
  _gradVx(),
  _gradVy(),
  //_gradVz(),
  socket_gravity("gravity"),
  socket_heating("Manchester"),
  socket_radiativeloss("RadiativeLossTerm"),
 // socket_zp("zplus"),
 // socket_zm("zminus"),
  socket_wavepressure("wavepressure"),
  socket_divBCellCenter("divBCellCenter")
{
  addConfigOptionsTo(this);
  _gravity = 0; 
  setParameter("gravity",&_gravity);
  setParameter("PevtsovHeating",&_PevtsovHeating);
  setParameter("PevtsovHeatingFactor",&_PevtsovHeatingFactor);
  _Manchester = 0;
  setParameter("Manchester",&_Manchester);
  setParameter("ManchesterHeatingAmplitude",&_ManchesterHeatingAmplitude);
  setParameter("ManchesterSigma",&_ManchesterSigma);
  setParameter("Qh4H_const", &_Qh4H_const);
  setParameter("Qlio_AR", &_Qlio_AR);
  setParameter("P0_W", &_P0_W);
  setParameter("alfven_pressure", &_alfven_pressure);
  setParameter("Qh2_activate", &_Qh2_activate);
  setParameter("Qh3_activate", &_Qh3_activate);
  setParameter("Qh4_activate", &_Qh4_activate);
  setParameter("Qh_lio_activate", &_Qh_lio_activate);
  setParameter("divQ",&_divQ);
  setParameter("divQConductivity",&_divQConductivity);
  setParameter("divQalphaCollisionless",&_divQalphaCollisionless);
  setParameter("Resistivity",&_Resistivity);
  setParameter("wavepressure", &_wave_pressure);  
  _RadiativeLossTerm = 0;
  setParameter("RadiativeLossTerm",&_RadiativeLossTerm);
  _2DdeCompE = 0;
  setParameter("2DdeCompEorNot", &_2DdeCompE);
  //_TriD = 1.0;
  //setParameter("TriD_factor", &_TriD);

//  setParameter("zplus",&_zplus);
//  setParameter("zminus",&_zminus); 
  //_projectionIDs = vector<CFuint>();
  //setParameter("ProjectionIDs",&_projectionIDs);

}

//////////////////////////////////////////////////////////////////////////////

void MHD2DConsACAHWSourceTerm::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFint >("gravity","Switch on gravity.");
  options.addConfigOption< CFint >("PevtsovHeating","Switch on PevtsovHeating term.");
  options.addConfigOption< CFreal >("PevtsovHeatingFactor","Scale PevtsovHeating term.");
  options.addConfigOption< CFint >("Manchester","Switch on Manchester heating term.");
  options.addConfigOption< CFreal >("ManchesterHeatingAmplitude","Scale Manchester volumetric heating amplitude.");
  options.addConfigOption< CFreal >("Qh4H_const", "Heating coefficient for the Qh4 heating function");
  options.addConfigOption< CFreal >("Qlio_AR", "Heating coefficient for the Qhlio heating function");
  options.addConfigOption< CFreal >("P0_W","amplitude for wave pressure");
  options.addConfigOption< CFint >("alfven_pressure", "flag for switching on alfven pressure");
  options.addConfigOption< CFint >("Qh2_activate", "whether Qh2 heating term should be activated");
  options.addConfigOption< CFint >("Qh3_activate", "whether Qh3 heating term should be activated");
  options.addConfigOption< CFint >("Qh4_activate", "whether Qh4 heating term should be activated");
  options.addConfigOption< CFint >("Qh_lio_activate", "whether Qh_lio heating term should be activated");
  options.addConfigOption< CFreal >("ManchesterSigma","Sigma value of Manchester heating term.");
  options.addConfigOption< CFint >("divQ","Switch on heat conduction.");
  options.addConfigOption< CFreal >("divQConductivity","Conductivity of heat conduction.");
  options.addConfigOption< CFreal >("divQalphaCollisionless","alpha value of heat conductivity in collisionless regime.");
  options.addConfigOption< CFreal >("Resistivity","eta value.");
  options.addConfigOption< CFint >("RadiativeLossTerm","Switch on optically thin approximation for radiation losses.");
  options.addConfigOption< CFint >("2DdeCompEorNot", "Calculate source term caused by decomposed energy equation.");
  //options.addConfigOption< CFreal >("TriD_factor", "Amount of fraction in Z- direction considered.");
  options.addConfigOption< CFint >("wavepressure", "save approximated wavepressure term.");
//  options.addConfigOption< CFint >("zplus", "z plus term from WKB");
//  options.addConfigOption< CFint >("zminus", "z minus term from WKB");
  
}

//////////////////////////////////////////////////////////////////////////////

MHD2DConsACAHWSourceTerm::~MHD2DConsACAHWSourceTerm()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DConsACAHWSourceTerm::setup()
{
  ComputeSourceTermFVMCC::setup();
  socket_gravity.getDataHandle().resize(socket_volumes.getDataHandle().size()*4);
  
  socket_heating.getDataHandle().resize(socket_volumes.getDataHandle().size());
  socket_radiativeloss.getDataHandle().resize(socket_volumes.getDataHandle().size());
  socket_wavepressure.getDataHandle().resize(socket_volumes.getDataHandle().size());
  socket_divBCellCenter.getDataHandle().resize(socket_volumes.getDataHandle().size());
  
  DataHandle<CFreal> normals = this->socket_normals.getDataHandle();
  _gradP.resize(PhysicalModelStack::getActive()->getDim(),0.);
  _gradBx.resize(PhysicalModelStack::getActive()->getDim(),0.);
  _gradBy.resize(PhysicalModelStack::getActive()->getDim(),0.);
  //_gradBz.resize(PhysicalModelStack::getActive()->getDim(),0.);
  _gradRho.resize(PhysicalModelStack::getActive()->getDim(),0.);
  _gradVx.resize(PhysicalModelStack::getActive()->getDim(),0.);
  _gradVy.resize(PhysicalModelStack::getActive()->getDim(),0.);
  //_gradVz.resize(PhysicalModelStack::getActive()->getDim(),0.);
}

//////////////////////////////////////////////////////////////////////////////

void MHD2DConsACAHWSourceTerm::computeSource(Framework::GeometricEntity *const element,
					 RealVector& source,
					 RealMatrix& jacobian)
{
  CFLog(DEBUG_MAX, "MHD2DConsACAHWSourceTerm::computeSource() => START\n");
  DataHandle<CFreal> normals = this->socket_normals.getDataHandle(); 
  DataHandle<CFreal> volumes = socket_volumes.getDataHandle();

  SafePtr<MHDTerm> model = PhysicalModelStack::getActive()->getImplementor()->
	                  getConvectiveTerm().d_castTo<MHDTerm>();

  // this source term is for MHD flows
  const vector<State*>* const states = element->getStates();
  const CFuint elementID = element->getID();
  
  // all elements in FVM should have only one state
  cf_assert(states->size() == 1);

  State *const currState = (*states)[0];
  
  const CFreal refSpeed = model->getRefSpeed();
  const CFreal refSpeedSq = refSpeed*refSpeed;

  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();
  const CFuint dim = PhysicalModelStack::getActive()->getDim(); 
  const CFuint startID = elementID*dim;

  for (CFuint i = 0; i < (nbEqs-1); ++i) {
    source[i] = 0.0;
  }

  const std::string correctionType = model->getCorrectionType();

  if (correctionType == "Mixed") {
    // mixed (hyperbolic and parabolic) correction
    const CFreal dissipCoeff = model->getDissipationCoefficient();
    const CFreal dissipCoeffSq = dissipCoeff*dissipCoeff;

    source[nbEqs-1] = -(refSpeedSq/dissipCoeffSq)*
      (*currState)[nbEqs-1]*
      volumes[elementID];
  }
  else {
    // hyperbolic correction
    source[nbEqs-1] = 0.0;
  }
  

 
  CFreal Bx =(*currState)[4]; 
  CFreal By =(*currState)[5];
  CFreal Bnorm = std::sqrt(std::pow(Bx,2) + std::pow(By,2));
  CFreal Vx = (*currState)[1];
  CFreal Vy =(*currState)[2];
	  //>> energy Source terms caused by decomposed energy equation
	  //CFuint deCompE = 1;
	  CFreal Q_deCompE = 0.0;
	  //if (deCompE == 1 && r_dimless>1.05){
	  //_deCompE = 1;
	  if (_2DdeCompE == 1){
		  const CFuint BXID = 4;
		  const CFuint BYID = 5;
		  const CFuint BZID = 6;
		  const CFuint gradBXID = elementID*nbEqs + BXID;
		  const CFuint gradBYID = elementID*nbEqs + BYID;
		  const CFuint gradBZID = elementID*nbEqs + BZID;
		  _gradBx[XX] = this->m_ux[gradBXID];
		  _gradBx[YY] = this->m_uy[gradBXID];
		  _gradBy[XX] = this->m_ux[gradBYID];
		  _gradBy[YY] = this->m_uy[gradBYID];
		  std::vector<CFreal> BdotLamdaB(2, 0.0);
		  std::vector<CFreal> VdotLamdaB(2, 0.0);
		  VdotLamdaB[0] = Vx*_gradBx[XX] + Vy*_gradBx[YY];
		  VdotLamdaB[1] = Vx*_gradBy[XX] + Vy*_gradBy[YY];
		  BdotLamdaB[0] = (Bx*_gradBx[XX] + By*_gradBx[YY]);
		  BdotLamdaB[1] = (Bx*_gradBy[XX] + By*_gradBy[YY]);
		  Q_deCompE = Vx*BdotLamdaB[0] + Vy*BdotLamdaB[1] -
			  (Bx*VdotLamdaB[0] + By*VdotLamdaB[1]);
		  source[7] += Q_deCompE*volumes[elementID];
	  }
	  //<< energy Source terms caused by decomposed energy equation

      source[8] = 0.0;
      // Eps modification 
  //    if (states->size() == 11) {
  //    source[9] = source_plus/vRef*volumes[elementID];
  //    source[10] = source_minus/vRef*volumes[elementID];
  //    socket_zp.getDataHandle()[elementID] = source_plus/vRef;
  //    socket_zm.getDataHandle()[elementID] = source_minus/vRef;
  //    }
  CFLog(DEBUG_MAX, "MHD2DConsACAHWSourceTerm::computeSource() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<Framework::BaseDataSocketSource> >
MHD2DConsACAHWSourceTerm::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result = ComputeSourceTermFVMCC::providesSockets();
  result.push_back(&socket_gravity);
  result.push_back(&socket_heating);
  result.push_back(&socket_radiativeloss);
  result.push_back(&socket_wavepressure);
  result.push_back(&socket_divBCellCenter);
  
  //  result.push_back(&socket_zp);
//  result.push_back(&socket_zm);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
