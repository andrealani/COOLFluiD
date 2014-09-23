#include "StructMech/StructMech.hh"
#include "StructMech3D.hh"
#include "StructMech3DDiffusiveDisp.hh"
#include "StructMech3DDiffusiveVarSet.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<StructMech3DDiffusiveDisp, DiffusiveVarSet, StructMechModule, 2> structMech3DDiffusiveDispProvider("StructMech3DDiffusiveDisp");

//////////////////////////////////////////////////////////////////////////////

void StructMech3DDiffusiveDisp::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< bool >("NonLinear","Option to turn on Geometrical NonLinearity if false will use Linear formulation.");
  options.addConfigOption< bool >("MeshMovement","Non uniform stiffness used for mesh movement.");
  options.addConfigOption< std::string >("MeshMovementMethod","Which kind of method to modify stiffness used for mesh movement.");

}

//////////////////////////////////////////////////////////////////////////////

StructMech3DDiffusiveDisp::StructMech3DDiffusiveDisp(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
  StructMech3DDiffusiveVarSet(name, model)
{
   addConfigOptionsTo(this);
  vector<std::string> names(3);
  names[0] = "u";
  names[1] = "v";
  names[2] = "w";
  setVarNames(names);

  _isNonLinear = false;
   setParameter("NonLinear",&_isNonLinear);

  _meshMovement = false;
   setParameter("MeshMovement",&_meshMovement);

  _meshMovementMethod = "VolumeBased";
   setParameter("MeshMovementMethod",&_meshMovementMethod);

}

//////////////////////////////////////////////////////////////////////////////

StructMech3DDiffusiveDisp::~StructMech3DDiffusiveDisp()
{
}

//////////////////////////////////////////////////////////////////////////////

void StructMech3DDiffusiveDisp::configure ( Config::ConfigArgs& args )
{
  StructMech3DDiffusiveVarSet::configure(args);
  bool isAnisotropic = getModel()->isAnisotropic();

  if (!isAnisotropic){
  // isotropic material
  const CFreal E  = getModel()->getYoung();
  const CFreal nu = getModel()->getPoisson();

  const CFreal mnu = 1. - nu;
  const CFreal D = E / ((1.+nu)*(1.-2.*nu));

  //Isotropic
  c11 = D * mnu;
  c12 = D * nu;
  c13 = D * nu;

  c21 = D * nu;
  c22 = D * mnu;
  c23 = D * nu;

  c31 = D * nu;
  c32 = D * nu;
  c33 = D * mnu;

  c44 = D * (1.0 - 2.0*nu)/2.0;
  c55 = D * (1.0 - 2.0*nu)/2.0;
  c66 = D * (1.0 - 2.0*nu)/2.0;

  //Set the other terms to zero
  c14 = 0.;
  c15 = 0.;
  c16 = 0.;
  c24 = 0.;
  c25 = 0.;
  c26 = 0.;
  c34 = 0.;
  c35 = 0.;
  c36 = 0.;
  c41 = 0.;
  c42 = 0.;
  c43 = 0.;
  c45 = 0.;
  c46 = 0.;
  c51 = 0.;
  c52 = 0.;
  c53 = 0.;
  c54 = 0.;
  c56 = 0.;
  c61 = 0.;
  c62 = 0.;
  c63 = 0.;
  c64 = 0.;
  c65 = 0.;

  /// Compute Lame Parameters from E and nu
  _lambda = E*nu/((1.+nu)*(1.-2*nu));
  _mu = E/(2*(1+nu));
  }
  else{

  ///@todo this should be done differently
  RealMatrix _c = *(getModel()->getStiffnessMatrix());

  cf_assert (_c.nbRows() == _c.nbCols());
  cf_assert (_c.nbRows() == 6);

  c11 = _c(0,0);
  c12 = _c(0,1);
  c13 = _c(0,2);
  c14 = _c(0,3);
  c15 = _c(0,4);
  c16 = _c(0,5);

  c21 = _c(1,0);
  c22 = _c(1,1);
  c23 = _c(1,2);
  c24 = _c(1,3);
  c25 = _c(1,4);
  c26 = _c(1,5);

  c31 = _c(2,0);
  c32 = _c(2,1);
  c33 = _c(2,2);
  c34 = _c(2,3);
  c35 = _c(2,4);
  c36 = _c(2,5);

  c41 = _c(3,0);
  c42 = _c(3,1);
  c43 = _c(3,2);
  c44 = _c(3,3);
  c45 = _c(3,4);
  c46 = _c(3,5);

  c51 = _c(4,0);
  c52 = _c(4,1);
  c53 = _c(4,2);
  c54 = _c(4,3);
  c55 = _c(4,4);
  c56 = _c(4,5);

  c61 = _c(5,0);
  c62 = _c(5,1);
  c63 = _c(5,2);
  c64 = _c(5,3);
  c65 = _c(5,4);
  c66 = _c(5,5);
  }

///Set the Stiffness Matrix
_stiffness.resize(6,6);
_stiffness(0,0) = c11;
_stiffness(0,1) = c12;
_stiffness(0,2) = c13;
_stiffness(0,3) = c14;
_stiffness(0,4) = c15;
_stiffness(0,5) = c16;
_stiffness(1,0) = c21;
_stiffness(1,1) = c22;
_stiffness(1,2) = c23;
_stiffness(1,3) = c24;
_stiffness(1,4) = c25;
_stiffness(1,5) = c26;
_stiffness(2,0) = c31;
_stiffness(2,1) = c32;
_stiffness(2,2) = c33;
_stiffness(2,3) = c34;
_stiffness(2,4) = c35;
_stiffness(2,5) = c36;
_stiffness(3,0) = c41;
_stiffness(3,1) = c42;
_stiffness(3,2) = c43;
_stiffness(3,3) = c44;
_stiffness(3,4) = c45;
_stiffness(3,5) = c46;
_stiffness(4,0) = c51;
_stiffness(4,1) = c52;
_stiffness(4,2) = c53;
_stiffness(4,3) = c54;
_stiffness(4,4) = c55;
_stiffness(4,5) = c56;
_stiffness(5,0) = c61;
_stiffness(5,1) = c62;
_stiffness(5,2) = c63;
_stiffness(5,3) = c64;
_stiffness(5,4) = c65;
_stiffness(5,5) = c66;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
