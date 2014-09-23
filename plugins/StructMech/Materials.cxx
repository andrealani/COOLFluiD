#include "Materials.hh"
#include "MaterialData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

Materials::Materials():
  _data(7),
  _temp()
{
addSteel();
addAluminium();
addTitanium();
addComposites();

_data.sortKeys();
}

//////////////////////////////////////////////////////////////////////////////

Materials::~Materials()
{
}

//////////////////////////////////////////////////////////////////////////////

void Materials::addSteel()
{
/// Data from "Mecanique des Materiaux",Massonnet & Cescoto pp366
/// Steel
_temp.name        = "Steel";
_temp.anisotropic = false;
_temp.C           = 0.;
_temp.young       = 205E9;
_temp.poisson     = 0.3;
_temp.density     = 7850.;
_temp.alpha       = 0.000012; //wikipedia, dilatation thermique
_data.insert(_temp.name,_temp);

_temp.name        = "Steel2";
_temp.anisotropic = false;
_temp.C           = 0.;
_temp.young       = 210E9;
_temp.poisson     = 0.28;
_temp.density     = 7700.;
_temp.alpha       = 0.000012; //wikipedia, dilatation thermique
_data.insert(_temp.name,_temp);

}

//////////////////////////////////////////////////////////////////////////////

void Materials::addAluminium()
{
/// Aluminium
_temp.name        = "Aluminium";
_temp.anisotropic = false;
_temp.C           = 0.;
_temp.young       = 65E9;
_temp.poisson     = 0.33;
_temp.density     = 2700.;
_temp.alpha       = 0.0000238; //wikipedia, dilatation thermique
_data.insert(_temp.name,_temp);

/// Aluminium
_temp.name        = "Al6061";
_temp.anisotropic = false;
_temp.C           = 0.;
_temp.young       = 70E9;
_temp.poisson     = 0.33;
_temp.density     = 2700.;
_temp.alpha       = 0.0000238; //wikipedia, dilatation thermique
_data.insert(_temp.name,_temp);

}

//////////////////////////////////////////////////////////////////////////////

void Materials::addTitanium()
{

/// Titane
_temp.name        = "Titanium";
_temp.anisotropic = false;
_temp.C           = 0.;
_temp.young       = 110E9;
_temp.poisson     = 0.25;
_temp.density     = 4500.;
_temp.alpha       = 0.0000105; //wikipedia, dilatation thermique
_data.insert(_temp.name,_temp);
}

//////////////////////////////////////////////////////////////////////////////

void Materials::addComposites()
{
/// Kevlar
_temp.name        = "Kevlar";
_temp.anisotropic = true;
RealMatrix C(6,6);
// Abaqus manual (taken from Kawabata and al. (1993))
// E and G values are in N/m2
CFreal E1 = 129600E6;
CFreal E2 = 2490E6;
CFreal E3 = 2490E6;

CFreal nu12 = 0.31;
CFreal nu23 = 0.0119;
CFreal nu13 = 0.62;

CFreal G12 = 2010E6;
CFreal G23 = 924E6;
CFreal G13 = 2010E6;

// First set all values to zero
C = 0;

// Only the Compliance matrix is known
// so invert it to get the stiffness matrix
RealMatrix S(3,3);
RealMatrix inverse(3,3);

// Set non-zero values
S(0,0)=1/E1;
S(0,1)=-nu12/E1;
S(0,2)=-nu13/E1;
S(1,0)=S(0,1);
S(1,1)=1/E2;
S(1,2)=-nu23/E2;
S(2,0)=S(0,2);
S(2,1)=S(1,2);
S(2,2)=1/E3;

///@todo change this, very bad (determinant is <<)
m_inverter3.invert(S,inverse);
for (CFuint i=0;i<3;++i){
  for (CFuint j=0;j<3;++j){
  C(i,j) = inverse(i,j);
  }
}

C(3,3)=G12;
C(4,4)=G23;
C(5,5)=G13;

_temp.C           = C;
_temp.young       = 0.;
_temp.poisson     = 0.;
_temp.density     = 1.0;
_data.insert(_temp.name,_temp);

/// Carbon/Epoxy (AS4/3501-6)
_temp.name        = "Carbon/Epoxy";
_temp.anisotropic = true;

// E and G values are in N/m2
E1 = 142E9;
E2 = 10.3E9;

nu12 = 0.27;
CFreal nu21 = 0.27;

G12 = 7.2E9;

nu13 = 0.;
CFreal nu31 = 0.;
nu23 = 0.;
CFreal nu32 = 0.;

E3 = 1.;
G13 = 1.;
G23 = 1.;

// First set all values to zero
C = 0;

// Only the Compliance matrix is known
// so invert it to get the stiffness matrix

// Set non-zero values
S(0,0)=1/E1;
S(0,1)=-nu21/E2;
S(0,2)=-nu31/E3;
S(1,0)=-nu21/E1;
S(1,1)=1/E2;
S(1,2)=-nu32/E3;
S(2,0)=-nu13/E1;
S(2,1)=-nu23/E2;
S(2,2)=1/E3;

///@todo change this, very bad (determinant is <<)
m_inverter3.invert(S,inverse);
for (CFuint i=0;i<3;++i){
  for (CFuint j=0;j<3;++j){
  C(i,j) = inverse(i,j);
  }
}

C(3,3)=G12;
C(4,4)=G23;
C(5,5)=G13;

_temp.C           = C;
_temp.young       = 0.;
_temp.poisson     = 0.;
_temp.density     = 1.0;
_data.insert(_temp.name,_temp);

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
