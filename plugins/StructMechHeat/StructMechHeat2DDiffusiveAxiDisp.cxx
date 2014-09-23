#include "StructMechHeat/StructMechHeat.hh"
#include "StructMechHeat2DDiffusiveAxiDisp.hh"
#include "StructMechHeat2DDiffusiveVarSet.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMechHeat {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<StructMechHeat2DDiffusiveAxiDisp, DiffusiveVarSet, StructMechHeatModule, 2> StructMechHeat2DDiffusiveAxiDispProvider("StructMechHeat2DDiffusiveAxiDisp");

//////////////////////////////////////////////////////////////////////////////

void StructMechHeat2DDiffusiveAxiDisp::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("PlaneStress","Option to turn on PlaneStress if false will use Plane Strain.");
  options.addConfigOption< bool >("NonLinear","Option to turn on Geometrical NonLinearity if false will use Linear formulation.");
  options.addConfigOption< bool >("MeshMovement","Non uniform stiffness used for mesh movement.");
  options.addConfigOption< std::string >("MeshMovementMethod","Which kind of method to modify stiffness used for mesh movement.");

}

//////////////////////////////////////////////////////////////////////////////

StructMechHeat2DDiffusiveAxiDisp::StructMechHeat2DDiffusiveAxiDisp(const std::string& name, Common::SafePtr<Framework::PhysicalModelImpl> model) :
  StructMechHeat2DDiffusiveVarSet(name, model)
{
   addConfigOptionsTo(this);
  vector<std::string> names(3);
  names[0] = "u";
  names[1] = "v";
  names[2] = "T";
  setVarNames(names);

  _isPlaneStress = false;
   setParameter("PlaneStress",&_isPlaneStress);

  _isNonLinear = false;
   setParameter("NonLinear",&_isNonLinear);

  _meshMovement = false;
   setParameter("MeshMovement",&_meshMovement);

  _meshMovementMethod = "VolumeBased";
   setParameter("MeshMovementMethod",&_meshMovementMethod);

}

//////////////////////////////////////////////////////////////////////////////

StructMechHeat2DDiffusiveAxiDisp::~StructMechHeat2DDiffusiveAxiDisp()
{
}

//////////////////////////////////////////////////////////////////////////////

void StructMechHeat2DDiffusiveAxiDisp::configure ( Config::ConfigArgs& args )
{

  CFAUTOTRACE;

  StructMechHeat2DDiffusiveVarSet::configure(args);

  ///Resize the Stiffness Matrix
  _stiffness.resize(4,4);
  _stiffness = 0.;
  // For an isotropic material,
  // Young:     Ex = Ey  = E
  // Poisson:   nuxy = nuyx = nu
  bool isAnisotropic = getModel()->isAnisotropic();
  if(!isAnisotropic){
    const CFreal E  = getModel()->getYoung();
    const CFreal nu = getModel()->getPoisson();
    const CFreal G = E /(2.*(1.+nu));

    if(_isPlaneStress) {
      cf_assert(false);
      // plane stress
    }
    else {
      // plane strain
      CFreal temp = E/((1.+nu)*(1.-2.*nu));
      _stiffness(0,0) = temp * (1.-nu);
      _stiffness(1,1) = temp * (1.-nu);
      _stiffness(2,2) = temp * (1.-nu);
      _stiffness(0,1) = temp * nu;
      _stiffness(0,2) = temp * nu;
      _stiffness(1,0) = temp * nu;
      _stiffness(1,2) = temp * nu;
      _stiffness(2,0) = temp * nu;
      _stiffness(2,1) = temp * nu;
      _stiffness(3,3) = G;

      /// Compute Lame Parameters from E and nu
      _lambda = E * nu;
      _lambda /= (1.+ nu)*(1.- 2.*nu);

      _mu = G;
    }
  }
  else{
    cf_assert(false);
  }

}

//////////////////////////////////////////////////////////////////////////////

RealMatrix& StructMechHeat2DDiffusiveAxiDisp::getStiffnessMat()
{
  return _stiffness;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMechHeat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
