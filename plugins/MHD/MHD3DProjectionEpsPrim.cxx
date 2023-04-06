#include "MHD/MHD.hh"

#include "MHD3DProjectionEpsPrim.hh"
#include "Common/CFLog.hh"
#include "Framework/State.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"
#include "Environment/ObjectProvider.hh"
#include "MathTools/MathFunctions.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MHD3DProjectionEpsPrim, ConvectiveVarSet, MHDModule, 1> 
mhd3DProjectionEpsPrimProvider("MHD3DProjectionEpsPrim");

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionEpsPrim::MHD3DProjectionEpsPrim(Common::SafePtr<Framework::BaseTerm> term) :
  MHD3DProjectionEpsVarSet(term)
{
  vector<std::string> names(11);
  names[0] = "rho";
  names[1] = "u";
  names[2] = "v";
  names[3] = "w";
  names[4] = "Bx";
  names[5] = "By";
  names[6] = "Bz";
  names[7] = "p";
  names[8] = "phi";
  names[9]  = "epsP";
  names[10] = "epsM";
  
  setVarNames(names);

  turnOnSourceTerm();
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionEpsPrim::~MHD3DProjectionEpsPrim()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionEpsPrim::setup()
{
  MHD3DProjectionEpsVarSet::setup();
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> MHD3DProjectionEpsPrim::getExtraVarNames() const
{
  vector<std::string> names(9);
  names[0] = "BxPotential";
  names[1] = "ByPotential";
  names[2] = "BzPotential";
  names[3] = "BxTotal";
  names[4] = "ByTotal";
  names[5] = "BzTotal";
  names[6] = "BTotal";
  names[7] = "rhoETotal";
  names[8] = "T";
    
  return names;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionEpsPrim::computeJacobians()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionEpsPrim::computeEigenValuesVectors(RealMatrix& rightEv,
						       RealMatrix& leftEv,
						       RealVector& eValues,
						       const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"MHD3DProjectionEpsPrim::computeEigenValuesVectors()");
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionEpsPrim::setDimensionalValuesPlusExtraValues(const State& state,
                                                              RealVector& result,
                                                              RealVector& extra)
{
  CFreal Bx0=0.0, By0=0.0, Bz0=0.0;
  const CFreal gammaMinus1 = getModel()->getGamma() - 1.;
  const CFreal B1dotB0 = state[4]*Bx0 + state[5]*By0 + state[6]*Bz0;
  const CFreal sqB0 = Bx0*Bx0 + By0*By0 + Bz0*Bz0;

  const CFreal rho = state[0];
  const CFreal u = state[1];
  const CFreal v = state[2];
  const CFreal w = state[3];
  const CFreal Bx1 = state[4];
  const CFreal By1 = state[5];
  const CFreal Bz1 = state[6];
  const CFreal V2 = u*u + v*v + w*w;
  const CFreal sqB1 = Bx1*Bx1 + By1*By1 + Bz1*Bz1;
  
  result[0] = state[0];
  result[1] = state[1];
  result[2] = state[2];
  result[3] = state[3];
  result[4] = state[4];
  result[5] = state[5];
  result[6] = state[6];
  result[7] = state[7];
  result[8] = state[8];
  result[9] = state[9];
  result[10] = state[10];
  
  extra.resize(9);
  
  extra[0] = Bx0;
  extra[1] = By0;
  extra[2] = Bz0;
  extra[3] = state[4] + extra[0];
  extra[4] = state[5] + extra[1];
  extra[5] = state[6] + extra[2];
  extra[6] = sqrt(extra[3]*extra[3] + extra[4]*extra[4] + extra[5]*extra[5]);
  extra[7] = (state[7]/gammaMinus1) + 0.5*(rho*V2+sqB1) + B1dotB0 + 0.5*sqB0;
  
  // AL: calculate temperature here from perfect gas
  extra[8] = 0.0 ; // 
}
      
//////////////////////////////////////////////////////////////////////////////

CFuint MHD3DProjectionEpsPrim::getBlockSeparator() const
{
  return PhysicalModelStack::getActive()->getNbEq();
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionEpsPrim::splitJacobian(RealMatrix& jacobPlus,
         RealMatrix& jacobMin,
         RealVector& eValues,
         const RealVector& normal)
{
  throw Common::NotImplementedException (FromHere(),"MHD3DProjectionEpsPrim::splitJacobian()");
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionEpsPrim::computePhysicalData(const State& state,
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

  data[MHDProjectionTerm::VX] = u;
  data[MHDProjectionTerm::VY] = v;
  data[MHDProjectionTerm::VZ] = w;
  data[MHDProjectionTerm::BX] = Bx;
  data[MHDProjectionTerm::BY] = By;
  data[MHDProjectionTerm::BZ] = Bz;
  data[MHDProjectionTerm::RHO] = rho;
  data[MHDProjectionTerm::V] = sqrt(V2);
  data[MHDProjectionTerm::B] = sqrt(B2);
  data[MHDProjectionTerm::P] = state[7];
  data[MHDProjectionTerm::PHI] = state[8];
  data[MHDProjectionTerm::A] = sqrt(getModel()->getGamma()*
                                    data[MHDProjectionTerm::P]/rho);
  
  data[MHDProjectionTerm::GAMMA] = getModel()->getGamma();
  data[MHDProjectionEpsTerm::EPSP] = state[9];
  data[MHDProjectionEpsTerm::EPSM] = state[10];
  
  const RealVector& node = state.getCoordinates();
  data[MHDTerm::XP] =  node[XX];
  data[MHDTerm::YP] =  node[YY];
  data[MHDTerm::ZP] =  node[ZZ];
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DProjectionEpsPrim::computeStateFromPhysicalData(const RealVector& data,
						  State& state)
{
  state[0] = data[MHDProjectionTerm::RHO];
  state[1] = data[MHDProjectionTerm::VX];
  state[2] = data[MHDProjectionTerm::VY];
  state[3] = data[MHDProjectionTerm::VZ];
  state[4] = data[MHDProjectionTerm::BX];
  state[5] = data[MHDProjectionTerm::BY];
  state[6] = data[MHDProjectionTerm::BZ];
  state[7] = data[MHDProjectionTerm::P];
  state[8] = data[MHDProjectionTerm::PHI];
  state[9] = data[MHDProjectionEpsTerm::EPSP];
  state[10] = data[MHDProjectionEpsTerm::EPSM];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

