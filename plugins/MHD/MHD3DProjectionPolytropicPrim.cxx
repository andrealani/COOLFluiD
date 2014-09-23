#include "MHD/MHD.hh"

#include "MHD3DProjectionPolytropicPrim.hh"
#include "Common/CFLog.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MHD3DProjectionPolytropicPrim, ConvectiveVarSet, MHDModule, 1> 
mhd3DProjectionPolytropicPrimProvider("MHD3DProjectionPolytropicPrim");

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionPolytropicPrim::MHD3DProjectionPolytropicPrim(Common::SafePtr<Framework::BaseTerm> term) :
  MHD3DProjectionPolytropicVarSet(term)
{
  vector<std::string> names(8);
  names[0] = "rho";
  names[1] = "u";
  names[2] = "v";
  names[3] = "w";
  names[4] = "Bx";
  names[5] = "By";
  names[6] = "Bz";
  names[7] = "phi";
  setVarNames(names);

  turnOnSourceTerm();
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionPolytropicPrim::~MHD3DProjectionPolytropicPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPolytropicPrim::setup()
{
  MHD3DProjectionPolytropicVarSet::setup();
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> MHD3DProjectionPolytropicPrim::getExtraVarNames() const
{
  vector<std::string> names(9);
  names[0] = "BxPotential";
  names[1] = "ByPotential";
  names[2] = "BzPotential";
  names[3] = "BxTotal";
  names[4] = "ByTotal";
  names[5] = "BzTotal";
  names[6] = "BTotal";
  names[7] = "p";
  names[8] = "divB";

  return names;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPolytropicPrim::computeJacobians()
{
  throw Common::NotImplementedException (FromHere(),"MHD3DProjectionPolytropicPrim::computeJacobians()");
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPolytropicPrim::setDimensionalValuesPlusExtraValues(const State& state,
                                                              RealVector& result,
                                                              RealVector& extra)
{
  const CFreal n = getModel()->getPolytropicIndex();
  const CFreal rho = state[0];

  result[0] = state[0];
  result[1] = state[1];
  result[2] = state[2];
  result[3] = state[3];
  result[4] = state[4];
  result[5] = state[5];
  result[6] = state[6];
  result[7] = state[7];

  extra.resize(9);

  const std::string potentialBType = getModel()->getPotentialBType();

  CFreal Bx0, By0, Bz0;
  RealVector stateCoordsCartesian(PhysicalModelStack::getActive()->getDim()),
        stateCoordsSpherical(PhysicalModelStack::getActive()->getDim()),
        BPFSSCartesian(PhysicalModelStack::getActive()->getDim());
  RealMatrix carSphTransMat(PhysicalModelStack::getActive()->getDim(),PhysicalModelStack::getActive()->getDim()),
        sphCarTransMat(PhysicalModelStack::getActive()->getDim(),PhysicalModelStack::getActive()->getDim());

  if (potentialBType == "Dipole") {
     RealVector&  BDipole = getMagneticDipole(state.getCoordinates()[XX],
                                              state.getCoordinates()[YY],
                                              state.getCoordinates()[ZZ]);

     Bx0 = BDipole[0];
     By0 = BDipole[1];
     Bz0 = BDipole[2];
  }
  if (potentialBType == "PFSS") {
     stateCoordsCartesian[0] = state.getCoordinates()[XX];
     stateCoordsCartesian[1] = state.getCoordinates()[YY];
     stateCoordsCartesian[2] = state.getCoordinates()[ZZ];
     setTransformationMatrices(stateCoordsCartesian,stateCoordsSpherical,carSphTransMat,sphCarTransMat);
     computePFSSMagneticField(stateCoordsSpherical,BPFSSCartesian,sphCarTransMat);

     // Boltzmann constant
     const CFreal k = 1.3806503e-23;

     // magnetic permeability at vacuum
     const CFreal mu0 = 4.e-7*MathTools::MathConsts::CFrealPi();

     // mass of proton
     const CFreal mp = 1.67262158e-27;

     // mass of electron
     const CFreal me = 9.10938188e-31;

     const CFreal nRef = getModel()->getNRef();
     const CFreal TRef = getModel()->getTRef();

     const CFreal rhoRef = nRef*(mp+me);
     const CFreal vRef = sqrt(2.0*k*TRef/mp);

     const CFreal BRef = sqrt(mu0*rhoRef*vRef*vRef)*1.e4;

     Bx0 = BPFSSCartesian[0]/BRef;
     By0 = BPFSSCartesian[1]/BRef;
     Bz0 = BPFSSCartesian[2]/BRef;
  }

  extra[0] = Bx0;
  extra[1] = By0;
  extra[2] = Bz0;
  extra[3] = state[4] + Bx0;
  extra[4] = state[5] + By0;
  extra[5] = state[6] + Bz0;
  extra[6] = sqrt(extra[3]*extra[3] + extra[4]*extra[4] + extra[5]*extra[5]);
  extra[7] = pow(rho,n)/n;

  std::string nsp = MeshDataStack::getActive()->getPrimaryNamespace();
  std::string datahandleName = nsp + "_divBNodal";

  DataHandle<CFreal> divBNodal = MeshDataStack::getActive()->getDataStorage()->getData<CFreal>(datahandleName);

  const CFuint stateID = state.getLocalID();

  extra[8] = divBNodal[stateID];
}

//////////////////////////////////////////////////////////////////////////////

CFuint MHD3DProjectionPolytropicPrim::getBlockSeparator() const
{
  return PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPolytropicPrim::computePhysicalData(const State& state,
					      RealVector& data)
{
  const CFreal rho = state[0];
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal w = state[3];
  const CFreal Bx = state[4];
  const CFreal By = state[5];
  const CFreal Bz = state[6];
  const CFreal V2 = u*u + v*v + w*w;
  const CFreal B2 = Bx*Bx + By*By + Bz*Bz;

  data[MHDProjectionPolytropicTerm::VX] = u;
  data[MHDProjectionPolytropicTerm::VY] = v;
  data[MHDProjectionPolytropicTerm::VZ] = w;
  data[MHDProjectionPolytropicTerm::BX] = Bx;
  data[MHDProjectionPolytropicTerm::BY] = By;
  data[MHDProjectionPolytropicTerm::BZ] = Bz;
  data[MHDProjectionPolytropicTerm::RHO] = rho;
  data[MHDProjectionPolytropicTerm::V] = sqrt(V2);
  data[MHDProjectionPolytropicTerm::B] = sqrt(B2);
  data[MHDProjectionPolytropicTerm::N] = getModel()->getPolytropicIndex();
  data[MHDProjectionPolytropicTerm::P] = pow(rho,data[MHDProjectionPolytropicTerm::N])/data[MHDProjectionPolytropicTerm::N];
  data[MHDProjectionPolytropicTerm::PHI] = state[7];
  data[MHDProjectionPolytropicTerm::A] = sqrt(getModel()->getGamma()*
                                    data[MHDProjectionPolytropicTerm::P]/rho);
  
  data[MHDProjectionPolytropicTerm::GAMMA] = getModel()->getGamma();
  
  const RealVector& node = state.getCoordinates();
  data[MHDProjectionPolytropicTerm::XP] =  node[XX];
  data[MHDProjectionPolytropicTerm::YP] =  node[YY];
  data[MHDProjectionPolytropicTerm::ZP] =  node[ZZ];
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionPolytropicPrim::computeStateFromPhysicalData(const RealVector& data,
						  State& state)
{
  state[0] = data[MHDProjectionPolytropicTerm::RHO];
  state[1] = data[MHDProjectionPolytropicTerm::VX];
  state[2] = data[MHDProjectionPolytropicTerm::VY];
  state[3] = data[MHDProjectionPolytropicTerm::VZ];
  state[4] = data[MHDProjectionPolytropicTerm::BX];
  state[5] = data[MHDProjectionPolytropicTerm::BY];
  state[6] = data[MHDProjectionPolytropicTerm::BZ];
  state[7] = data[MHDProjectionPolytropicTerm::PHI];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

