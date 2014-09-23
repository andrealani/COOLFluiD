#include "StructMechHeat/StructMechHeat.hh"
#include "Framework/State.hh"
#include "StructMechHeat2D.hh"
#include "Environment/ObjectProvider.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace StructMechHeat {


//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<StructMechHeat2D,
               PhysicalModelImpl,
	       StructMechHeatModule,
               1>
StructMechHeat2DProvider("StructMechHeat2D");

//////////////////////////////////////////////////////////////////////////////

void StructMechHeat2D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Thickness","Thickness.");
  options.addConfigOption< bool >("IsAxisymmetric","IsAxisymmetric flag.");
  options.addConfigOption< std::string >("AxisymmetryAxis","Definition of the axis of revolution.");
}

//////////////////////////////////////////////////////////////////////////////

StructMechHeat2D::StructMechHeat2D(const std::string& name)
  : StructMechHeatPhysicalModel(name)
{
   addConfigOptionsTo(this);

   m_thickness = 1.0;
   setParameter("Thickness",&m_thickness);

   m_isAxisymmetric = false;
   setParameter("IsAxisymmetric",&m_isAxisymmetric);

   m_axisymmetryAxis = "X";
   setParameter("AxisymmetryAxis",&m_axisymmetryAxis);
}

//////////////////////////////////////////////////////////////////////////////

StructMechHeat2D::~StructMechHeat2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void StructMechHeat2D::configure ( Config::ConfigArgs& args )
{
  StructMechHeatPhysicalModel::configure(args);

  /// @note Should be here, but is not yet being used.
  _jacobians = std::vector<RealMatrix>(getDimension());
  for (CFuint i = 0; i < getDimension(); ++i) {
    _jacobians[i] = RealMatrix(getNbEquations(), getNbEquations());
    _jacobians[i] = 0.0;
  }

  if((m_axisymmetryAxis == "X") || (m_axisymmetryAxis == "x")){
    m_radiusID = YY;
  }
  else{
    cf_assert((m_axisymmetryAxis == "Y") || (m_axisymmetryAxis == "y"));
    m_radiusID = XX;
  }

  if(m_isAxisymmetric && (m_thickness != 1.) )
  {
    CFout << "WARNING: You have given a thickness in a Axisymmetric computation!!\n";
    CFout << "Thickness value will be ignored\n";
  }

  cf_assert(m_thickness > 0.);
  m_nbDim = getDimension();
  m_nbEqs = getNbEquations();

}

//////////////////////////////////////////////////////////////////////////////

CFuint StructMechHeat2D::getDimension() const
{
  return 2;
}

//////////////////////////////////////////////////////////////////////////////

CFuint StructMechHeat2D::getNbEquations() const
{
  return 3;
}

//////////////////////////////////////////////////////////////////////////////

bool StructMechHeat2D::validate(const State& state) const
{
  return state[2] >= 0.0;
}

//////////////////////////////////////////////////////////////////////////////

CFreal StructMechHeat2D::getThickness(const RealVector& coord)
{
  if(m_isAxisymmetric){
    m_thickness = 2. * MathTools::MathConsts::CFrealPi() * std::fabs(coord[m_radiusID]);
  }

  return m_thickness;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace StructMechHeat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
