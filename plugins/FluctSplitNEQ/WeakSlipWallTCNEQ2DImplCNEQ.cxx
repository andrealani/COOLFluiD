#include "FluctSplitNEQ/FluctSplitNEQ.hh"
#include "FluctSplitNEQ/WeakSlipWallTCNEQ2DImplCNEQ.hh"
#include "FluctSplit/InwardNormalsData.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::FluctSplit;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplitNEQ {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<WeakSlipWallTCNEQ2DImplCNEQ,
		      FluctuationSplitData, 
		      FluctSplitNEQModule> 
weakSlipWallTCNEQ2DImplCNEQProvider("WeakSlipWallTCNEQ2DImplCNEQ");

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallTCNEQ2DImplCNEQ::WeakSlipWallTCNEQ2DImplCNEQ(const std::string& name) :
   WeakSlipWall2DImpl(name),
   _varSet(CFNULL),
   _library(CFNULL),
   _ys(),
   _Ms()
{
}

//////////////////////////////////////////////////////////////////////////////

WeakSlipWallTCNEQ2DImplCNEQ::~WeakSlipWallTCNEQ2DImplCNEQ()
{
}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallTCNEQ2DImplCNEQ::setup()
{
  WeakSlipWall2DImpl::setup();
  
  _varSet = getMethodData().getUpdateVar().d_castTo<TCNEQ2DVarSet>();
  
  // set up the physical data
  _varSet->getModel()->resizePhysicalData(_physicalData);
  
  _library = PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  cf_assert(_library.isNotNull());
  
  const CFuint nbSpecies = _library->getNbSpecies();
  _ys.resize(nbSpecies);
  _Ms.resize(nbSpecies);
  
  // check the units
  _library->getMolarMasses(_Ms);
}

//////////////////////////////////////////////////////////////////////////////

void WeakSlipWallTCNEQ2DImplCNEQ::computeNormalFluxAndJacob
(const State& state,
 const RealVector& normal,
 RealVector& flux,
 RealMatrix& fluxJacob)
{
  State& ss = *(const_cast<State*>(&state));
  
  // check tthis !!! before it was not correct ...
  _varSet->setExtraData(true);
  _varSet->computePhysicalData(ss, _physicalData); 
  _varSet->setExtraData(false);
  
  const CFreal nx = normal[XX];
  const CFreal ny = normal[YY];
  const CFreal rho = _physicalData[EulerTerm::RHO];
  const CFreal u = _physicalData[EulerTerm::VX];
  const CFreal v = _physicalData[EulerTerm::VY];
  const CFreal q = 0.5*(u*u + v*v);		
  const CFreal un = u*nx + v*ny;
  const CFreal H = _physicalData[EulerTerm::H];
  const CFreal T = _physicalData[EulerTerm::T];
  SafePtr<PhysicalChemicalLibrary::ExtraData> eData = _library->getExtraData();
  
  // we assume to have [rhoi v T ] variables
  const CFuint nbSpecies = _varSet->getModel()->getNbScalarVars(0);
  CFreal ovM = 0.0;
  for (CFuint is = 0; is < nbSpecies; ++is) {
    flux[is] = state[is]*un;
    ovM += state[is]/_Ms[is];
  }
  
  const CFuint uID = nbSpecies; 
  const CFuint vID = nbSpecies+1; 
  const CFuint eID = nbSpecies+2;   

  const CFreal R = _library->getRgas();
  const CFreal rt = R*T;
  const CFreal ovRho = 1./rho;
  const CFreal pOvRho = _physicalData[EulerTerm::P]*ovRho;
  const CFreal hT = eData->dEdT + pOvRho/T; // check this !!!!

  flux[uID] = un*rho*u;
  flux[vID] = un*rho*v;
  flux[eID] = un*rho*H;

  // the analytical jacobian of the normal fluxes dFn/dP
  for (CFuint is = 0; is < nbSpecies; ++is) {
    fluxJacob(is,is) = un;
    fluxJacob(is,uID) = state[is]*nx;
    fluxJacob(is,vID) = state[is]*ny;
    
    fluxJacob(uID,is) = u*un;
    fluxJacob(vID,is) = v*un;
    
    // rho*H*un = un*(rho*E + p)
    const CFreal dRhoEdRhoi = eData->dRhoEdRhoi[is] + q;
    fluxJacob(eID,is) = un*dRhoEdRhoi + un*(rt/_Ms[is]);    
  }
  
  fluxJacob(uID,uID) = rho*(2.*u*nx + v*ny);
  fluxJacob(uID,vID) = rho*u*ny;
  
  fluxJacob(vID,uID) = rho*v*nx;
  fluxJacob(vID,vID) = rho*(u*nx + 2.*v*ny);
  
  fluxJacob(eID,uID) = rho*(H*nx + un*u);
  fluxJacob(eID,vID) = rho*(H*ny + un*v);
  fluxJacob(eID,eID) = rho*hT*un;

}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluctSplitNEQ



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
