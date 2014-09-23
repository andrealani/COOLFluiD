#include "Heat/Heat.hh"
#include "Framework/State.hh"
#include "Heat2D.hh"
#include "Environment/ObjectProvider.hh"
#include "MathTools/MathConsts.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace Heat {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Heat2D, PhysicalModelImpl, HeatModule, 1> heat2DProvider("Heat2D");

//////////////////////////////////////////////////////////////////////////////

void Heat2D::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("Thickness","Thickness.");
  options.addConfigOption< bool >("IsAxisymmetric","IsAxisymmetric flag.");
  options.addConfigOption< std::string >("AxisymmetryAxis","Definition of the axis of revolution.");
}

//////////////////////////////////////////////////////////////////////////////

Heat2D::Heat2D(const std::string& name)
  : HeatPhysicalModel(name)
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

Heat2D::~Heat2D()
{
}

//////////////////////////////////////////////////////////////////////////////

void Heat2D::configure ( Config::ConfigArgs& args )
{
  HeatPhysicalModel::configure(args);

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

  m_variables.resize(m_nbDim + m_nbEqs);

}

//////////////////////////////////////////////////////////////////////////////

CFuint Heat2D::getDimension() const
{
  return 2;
}

//////////////////////////////////////////////////////////////////////////////

CFuint Heat2D::getNbEquations() const
{
  return 1;
}

//////////////////////////////////////////////////////////////////////////////

bool Heat2D::validate(const State& state) const
{
  return state[0] >= 0.0;
}

//////////////////////////////////////////////////////////////////////////////

CFreal Heat2D::getThickness(const RealVector& coord)
{

  if(m_isAxisymmetric){
    m_thickness = 2. * MathTools::MathConsts::CFrealPi() * std::fabs(coord[m_radiusID] + MathTools::MathConsts::CFrealEps());
  }

  return m_thickness;

}

//////////////////////////////////////////////////////////////////////////////

    } // namespace Heat

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
