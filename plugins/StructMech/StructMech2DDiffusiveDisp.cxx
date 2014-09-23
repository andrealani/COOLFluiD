#include "StructMech/StructMech.hh"
#include "StructMech2DDiffusiveDisp.hh"
#include "StructMech2DDiffusiveVarSet.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMech {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<StructMech2DDiffusiveDisp, DiffusiveVarSet, StructMechModule, 2> structMech2DDiffusiveDispProvider("StructMech2DDiffusiveDisp");

//////////////////////////////////////////////////////////////////////////////

void StructMech2DDiffusiveDisp::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("PlaneStress","Option to turn on PlaneStress if false will use Plane Strain.");
  options.addConfigOption< bool >("NonLinear","Option to turn on Geometrical NonLinearity if false will use Linear formulation.");
  options.addConfigOption< bool >("MeshMovement","Non uniform stiffness used for mesh movement.");
  options.addConfigOption< std::string >("MeshMovementMethod","Which kind of method to modify stiffness used for mesh movement.");

}

//////////////////////////////////////////////////////////////////////////////

StructMech2DDiffusiveDisp::StructMech2DDiffusiveDisp(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
  StructMech2DDiffusiveVarSet(name, model)
{
   addConfigOptionsTo(this);
  vector<std::string> names(2);
  names[0] = "u";
  names[1] = "v";
  setVarNames(names);

  _isPlaneStress = true;
   setParameter("PlaneStress",&_isPlaneStress);

  _isNonLinear = false;
   setParameter("NonLinear",&_isNonLinear);

  _meshMovement = false;
   setParameter("MeshMovement",&_meshMovement);

  _meshMovementMethod = "VolumeBased";
   setParameter("MeshMovementMethod",&_meshMovementMethod);

}

//////////////////////////////////////////////////////////////////////////////

StructMech2DDiffusiveDisp::~StructMech2DDiffusiveDisp()
{
}

//////////////////////////////////////////////////////////////////////////////

void StructMech2DDiffusiveDisp::configure ( Config::ConfigArgs& args )
{

  CFAUTOTRACE;

  StructMech2DDiffusiveVarSet::configure(args);

  // For an isotropic material,
  // Young:     Ex = Ey  = E
  // Poisson:   nuxy = nuyx = nu
  bool isAnisotropic = getModel()->isAnisotropic();
  if(!isAnisotropic){
    const CFreal E  = getModel()->getYoung();
    const CFreal nu = getModel()->getPoisson();
    const CFreal G = E /(2.*(1.+nu));

    _c61 = 0.;
    _c62 = 0.;
    _c16 = 0.;
    _c26 = 0.;

    if(_isPlaneStress) {
      // plane stress
      _c11 = E/(1.-nu*nu);
      _c22 = _c11;
      _c12 = nu * _c11;
      _c21 = _c12;
      _c66 = G;
      /// Compute Lame Parameters from E and nu
      _lambda = _c11 * nu;
      _mu = G;
    }
    else {
      // plane strain
      _c11 = E*(1.-nu)/(1.-nu-2.*nu*nu);
      _c22 = E*(1.-nu*nu)/((1.-nu-2.*nu*nu)*(1.+nu));
      _c12 = E*nu/(1.-nu-2.*nu*nu);
      _c21 = _c12;
      _c66 = G;

      /// Compute Lame Parameters from E and nu
      _lambda = E * nu;
      _lambda /= (1.+ nu)*(1.- 2.*nu);

      _mu = G;
    }
  }
  else{
    /// Anisotropic Material
    /// @todo this should be done differently
    RealMatrix _c = *(getModel()->getStiffnessMatrix());

    cf_assert (_c.nbRows() == _c.nbCols());
    cf_assert (_c.nbRows() == 6);

    _c11 = _c(0,0);
    _c12 = _c(0,1);
    _c16 = _c(0,3);

    _c21 = _c(1,0);
    _c22 = _c(1,1);
    _c26 = _c(1,3);

    _c61 = _c(3,0);
    _c62 = _c(3,1);
    _c66 = _c(3,3);
  }

///Set the Stiffness Matrix
_stiffness.resize(3,3);
_stiffness(0,0) = _c11;
_stiffness(0,1) = _c12;
_stiffness(0,2) = _c16;
_stiffness(1,0) = _c21;
_stiffness(1,1) = _c22;
_stiffness(1,2) = _c26;
_stiffness(2,0) = _c61;
_stiffness(2,1) = _c62;
_stiffness(2,2) = _c66;

}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& StructMech2DDiffusiveDisp::getStiffnessMat()
{
  return _stiffness;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMech

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
