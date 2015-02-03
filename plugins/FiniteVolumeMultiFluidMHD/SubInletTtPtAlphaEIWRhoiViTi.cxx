#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/SubInletTtPtAlphaEIWRhoiViTi.hh"
#include "Framework/MethodCommandProvider.hh"
#include "MultiFluidMHD/DiffMFMHD2DVarSet.hh"
#include "MultiFluidMHD/MultiFluidMHDVarSet.hh"
#include "MultiFluidMHD/EulerMFMHDTerm.hh"
#include "Maxwell/Maxwell2DProjectionVarSet.hh"
#include "Framework/PhysicalConsts.hh"
// #include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::MultiFluidMHD;
using namespace COOLFluiD::Physics::Maxwell;
using namespace COOLFluiD::Config;


//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<SubInletTtPtAlphaEIWRhoiViTi, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule> subInletTtPtAlphaEIWRhoiViTiFVMCCProvider("SubInletTtPtAlphaEIWRhoiViTiFVMCC");


//////////////////////////////////////////////////////////////////////////////

void SubInletTtPtAlphaEIWRhoiViTi::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("Ttot","total temperature");
   options.addConfigOption< CFreal >("Ptot","total pressure");
   options.addConfigOption< CFreal >("alpha","alpha");
}

//////////////////////////////////////////////////////////////////////////////

SubInletTtPtAlphaEIWRhoiViTi::SubInletTtPtAlphaEIWRhoiViTi(const std::string& name) :
  FVMCC_BC(name),
  _varSet(CFNULL),
  _dataInnerState(),
  _dataGhostState()
{
   addConfigOptionsTo(this);
  _tTotal = 0.0;
   setParameter("Ttot",&_tTotal);

   _pTotal = 0.0;
   setParameter("Ptot",&_pTotal);

   _alpha = 0.0;
   setParameter("alpha",&_alpha);
}

//////////////////////////////////////////////////////////////////////////////

SubInletTtPtAlphaEIWRhoiViTi::~SubInletTtPtAlphaEIWRhoiViTi()
{
}

//////////////////////////////////////////////////////////////////////////////

void SubInletTtPtAlphaEIWRhoiViTi::setGhostState(GeometricEntity *const face)
{
    const CFuint nbSpecies = _varSet->getModel()->getNbScalarVars(0);

    State *const innerState = face->getState(0);
    State *const ghostState = face->getState(1);

    // set the physical data starting from the inner state
    _varSet->computePhysicalData(*innerState, _dataInnerState);


    ///Maxwell Equations Perfectly Conducting Wall Condition
    const CFuint faceID = face->getID();
    const CFuint startID = faceID*PhysicalModelStack::getActive()->getDim();

    DataHandle<CFreal> normals = socket_normals.getDataHandle();
    CFreal nx = normals[startID];
    CFreal ny = normals[startID + 1];
    CFreal nz = 0;
    const CFreal invFaceLength = 1./sqrt(nx*nx + ny*ny);
    nx *= invFaceLength;
    ny *= invFaceLength;

    cf_assert(_varSet.isNotNull());

    const CFreal bn = (_dataInnerState)[0]*nx + (_dataInnerState)[1]*ny;
    const CFreal en = (_dataInnerState)[3]*nx + (_dataInnerState)[4]*ny;
  //  const CFreal chi = _varSet->getModel()->getDivECleaningConst();

    (_dataGhostState)[0] = (_dataInnerState)[0] /*+ 2*bn*nx*/;	//Bx
    (_dataGhostState)[1] = (_dataInnerState)[1] /*+ 2*bn*ny*/;	//By
    (_dataGhostState)[2] = (_dataInnerState)[2] /*+ 2*bn*nz*/;	//Bz
    (_dataGhostState)[3] = (_dataInnerState)[3] /*- 2*en*nx*/;	//Ex
    (_dataGhostState)[4] = (_dataInnerState)[4] /*- 2*en*ny*/;	//Ey
    (_dataGhostState)[5] = (_dataInnerState)[5] /*- 2*en*nz*/;	//Ez
    (_dataGhostState)[6] = (_dataInnerState)[6];			//Psi
    (_dataGhostState)[7] = (_dataInnerState)[7];			//Phi

  ///MultiFluidMHD Subsonic Outlet Condition imposing pressure in 2D
   // here a fix is needed in order to have always ghostT > 0
   // if ghostT < 0  then the inner value is set
    const CFuint endEM = 8;
    const CFuint firstSpecies = _varSet->getModel()->getFirstScalarVar(0);
    const CFuint firstVelocity = _varSet->getModel()->getFirstScalarVar(1);
    const CFuint firstTemperature = _varSet->getModel()->getFirstScalarVar(2);
    const CFreal gamma = _varSet->getModel()->getGamma();
    const CFreal gammaDivGammaMinus1 = gamma/(gamma - 1.0);

    for (CFuint ie = 0; ie < nbSpecies; ++ie) {
        const CFreal m_air = _varSet->getModel()->getMolecularMass1();
        const CFreal K_B = PhysicalConsts::Boltzmann();
        const CFreal Rgas = K_B/m_air;				// ions gas constant

        const CFreal uxi = (_dataInnerState)[firstVelocity + 2*ie];
        const CFreal uyi = (_dataInnerState)[firstVelocity + 2*ie + 1];
        const CFreal Ui = std::sqrt(uxi*uxi + uyi*uyi);
        const CFreal rhoi = (_dataInnerState)[endEM]*(_dataInnerState)[firstSpecies + ie];
        const CFreal Ti = (_dataInnerState)[firstTemperature + 4*ie];
        const CFreal pi = rhoi*Rgas*Ti;
        const CFreal ai = (_dataInnerState)[firstTemperature + 4*ie + 2];
        const CFreal Mi = Ui/ai;
        const CFreal coeff = 1 + 0.5*(gamma - 1)*Mi*Mi;
        const CFreal Ttoti = Ti*coeff;
        //const CFreal coeff2 = (1 + gamma*Mi*Mi);
        const CFreal Ptoti = pi + 0.5*rhoi*Ui*Ui;
        const CFreal tgAlphai = uyi/uxi;

        const CFreal Ttotg = 2*_tTotal - Ttoti;
        const CFreal Ptotg = 2*_pTotal - Ptoti;
        const CFreal tgAlphag = 2*tan(_alpha) - tgAlphai;
        const CFreal Mg = Mi;
        const CFreal Tg = Ttotg/(1 + 0.5*(gamma - 1)*Mg*Mg);
        const CFreal pg = Ptotg/(1 + gamma*Mg*Mg);
        const CFreal rhog = pg/(Rgas*Tg);
        const CFreal uxg = Mg*std::sqrt(gamma*Rgas*Tg/(1 + tgAlphag*tgAlphag));
        const CFreal uyg = tgAlphag*uxg;
        const CFreal UgUg = uxg*uxg + uyg*uyg;
        const CFreal Hg = (gamma/(gamma - 1)*pg + 0.5*rhog*UgUg)/rhog;
        const CFreal Ag = sqrt(gamma*pg/rhog);


	//cout <<"m_air = " << m_air <<"\n";
        //cout <<"K_B = " << K_B <<"\n";
        //cout <<"Rgas = " << Rgas <<"\n";
        //cout <<"Ttotg " << Ttotg <<"\n";
        //cout <<"CFreal Ptotg = " << Ptotg <<"\n";
        //cout <<"tgAlphag = " << tgAlphag <<"\n";
        //cout <<"Mg = " << Mg <<"\n";
        //cout <<"Tg = " << Tg <<"\n";
        //cout <<"pg = " << pg <<"\n";
        //cout <<"rhog = " << rhog <<"\n";
        //cout <<"uxg = " << uxg <<"\n";
        //cout <<"uyg = " << uyg <<"\n";
        //cout <<"UgUg = " << UgUg <<"\n";
        //cout <<"Hg = " << Hg <<"\n";
        //cout <<"Ag = " << Ag <<"\n";

        (_dataGhostState)[endEM] = rhog;
        (_dataGhostState)[firstSpecies + ie] = 1; //Only one species;
        (_dataGhostState)[firstVelocity + 2*ie] = uxg;
        (_dataGhostState)[firstVelocity + 2*ie + 1] = uyg;
        (_dataGhostState)[firstTemperature + 4*ie] = Tg;
        (_dataGhostState)[firstTemperature + 4*ie + 1] = pg;
        (_dataGhostState)[firstTemperature + 4*ie + 2] = Ag;
        (_dataGhostState)[firstTemperature + 4*ie + 3] = Hg;
    }
    _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);


//    (_dataGhostState)[endEM] = (_dataInnerState)[endEM]; 							//RHO

//    for (CFuint ie = 0; ie < nbSpecies; ++ie) {
//      const CFreal Pi = _Pi[ie];
//  //     cout << "SubOutletPEIWRhoiViTi::setGhostState => Pi = " << Pi << endl;
//      (_dataGhostState)[firstSpecies + ie] = (_dataInnerState)[firstSpecies + ie];			//yi
//      (_dataGhostState)[firstVelocity + 2*ie] = (_dataInnerState)[firstVelocity + 2*ie];			//Vxi
//      (_dataGhostState)[firstVelocity + 2*ie + 1] = (_dataInnerState)[firstVelocity + 2*ie + 1];		//Vyi
//      (_dataGhostState)[firstTemperature + 4*ie] = (_dataInnerState)[firstTemperature + 4*ie];		//Ti
//      (_dataGhostState)[firstTemperature + 4*ie + 1] = 2*Pi -
//                              (_dataInnerState)[firstTemperature + 4*ie + 1];	//Pi

//      const CFreal Vi2 = (_dataGhostState)[firstVelocity + 2*ie]*(_dataGhostState)[firstVelocity + 2*ie] +
//              (_dataGhostState)[firstVelocity + 2*ie + 1]*(_dataGhostState)[firstVelocity + 2*ie + 1];
//      const CFreal rhoi =(_dataGhostState)[endEM]*(_dataGhostState)[firstSpecies + ie];

//      (_dataGhostState)[firstTemperature + 4*ie + 2] = sqrt(gamma*(_dataGhostState)[firstTemperature + 4*ie + 1]/rhoi);		//ai
//      (_dataGhostState)[firstTemperature + 4*ie + 3] = (0.5*rhoi*Vi2 +
//                              gammaDivGammaMinus1*(_dataGhostState)[firstTemperature + 4*ie + 1])/rhoi;	//Hi
//    }

//    _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);

//  // set the physical data starting from the inner state
//  _varSet->computePhysicalData(*innerState, _dataInnerState);

//  const CFreal gamma = _varSet->getModel()->getGamma();
//  const CFreal gammaMinus1 = gamma - 1.;
//  const CFreal gammaDivGammaMinus1 = gamma/gammaMinus1;
//  const CFreal R = _varSet->getModel()->getR();
//  const CFreal u = _dataInnerState[EulerTerm::VX];
//  const CFreal v = _dataInnerState[EulerTerm::VY];
//  const CFreal uSqvSq = u*u + v*v;
//  const CFreal pInnerState = _varSet->getModel()->getPressureFromState(_dataInnerState[EulerTerm::P]);
//  const CFreal tInnerState = pInnerState / (R*_dataInnerState[EulerTerm::RHO]);
//  const CFreal machInner = sqrt(uSqvSq/(gamma*R*tInnerState));
//  const CFreal coefficient = 1.0 + 0.5*gammaMinus1*machInner*machInner;
//  const CFreal tTotalInner = tInnerState*coefficient;
//  const CFreal coeffPow = pow(coefficient, gamma/gammaMinus1);
//  const CFreal pTotalInner = pInnerState*coeffPow;
//  const CFreal tgAlphaInner = v/u;
  
//  // ghost state quantities
//  const CFreal tTotalGhost = 2.0*_tTotal - tTotalInner;
//  const CFreal pTotalGhost = 2.0*_pTotal - pTotalInner;
//  const CFreal tgAlphaGhost = 2.0*tan(_alpha) - tgAlphaInner;
//  const CFreal machGhost = machInner;
//  const CFreal tGhost = tTotalGhost / coefficient;
  
//  // set the physical data for the ghost state
//  const CFreal pressure = pTotalGhost/coeffPow;
//  _dataGhostState[EulerTerm::P] = pressure - _varSet->getModel()->getPressInf();
//  _dataGhostState[EulerTerm::RHO] = pressure/(R*tGhost);
//  _dataGhostState[EulerTerm::VX] = machGhost*sqrt(gamma*R*tGhost/(1.0 + tgAlphaGhost*tgAlphaGhost));
//  _dataGhostState[EulerTerm::VY] = tgAlphaGhost*_dataGhostState[EulerTerm::VX];
//  _dataGhostState[EulerTerm::V] = sqrt(_dataGhostState[EulerTerm::VX]*
//				       _dataGhostState[EulerTerm::VX] +
//				       _dataGhostState[EulerTerm::VY]*
//				       _dataGhostState[EulerTerm::VY]);
//  _dataGhostState[EulerTerm::H] = (gammaDivGammaMinus1*pressure + 0.5*_dataGhostState[EulerTerm::RHO]*
//				   _dataGhostState[EulerTerm::V]*
//				   _dataGhostState[EulerTerm::V])/_dataGhostState[EulerTerm::RHO];
//  _dataGhostState[EulerTerm::A] = sqrt(gamma*pressure/_dataGhostState[EulerTerm::RHO]);
//  _dataGhostState[EulerTerm::T] = tGhost;
  
//  // set the ghost state starting from the physical data
//  _varSet->computeStateFromPhysicalData(_dataGhostState, *ghostState);
}

//////////////////////////////////////////////////////////////////////////////

void SubInletTtPtAlphaEIWRhoiViTi::setup()
{
  FVMCC_BC::setup();
	
  _varSet = getMethodData().getUpdateVar().d_castTo<MultiFluidMHDVarSet<Maxwell2DProjectionVarSet> >();
  _varSet->getModel()->resizePhysicalData(_dataInnerState);
  _varSet->getModel()->resizePhysicalData(_dataGhostState);

//  _pTotal/=_varSet->getModel()->getPressRef();
//  _tTotal/=_varSet->getModel()->getTempRef();

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
