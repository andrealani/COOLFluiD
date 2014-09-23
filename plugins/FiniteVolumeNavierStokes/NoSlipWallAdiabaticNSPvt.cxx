#include "FiniteVolumeNavierStokes/FiniteVolumeNavierStokes.hh"
#include "FiniteVolumeNavierStokes/NoSlipWallAdiabaticNSPvt.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"
#include "Framework/PhysicalModel.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NoSlipWallAdiabaticNSPvt, CellCenterFVMData, FiniteVolumeNavierStokesModule>
noSlipWallAdiabaticNSPvtFVMCCProvider("NoSlipWallAdiabaticNSPvtFVMCC");

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallAdiabaticNSPvt::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< std::vector<std::string> >
    ("WallVelocity","Definition of the functions defining the wall velocity.");
}

//////////////////////////////////////////////////////////////////////////////

NoSlipWallAdiabaticNSPvt::NoSlipWallAdiabaticNSPvt(const std::string& name) :
  FVMCC_BC(name)
{
  addConfigOptionsTo(this);

  m_functions = std::vector<std::string>();
  setParameter("WallVelocity",&m_functions);
}

//////////////////////////////////////////////////////////////////////////////

NoSlipWallAdiabaticNSPvt::~NoSlipWallAdiabaticNSPvt()
{
}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallAdiabaticNSPvt::configure ( Config::ConfigArgs& args )
{
  FVMCC_BC::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallAdiabaticNSPvt::setupWallVelocity()
{
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  // if the functions are empty them the velocity is set to zero
  if (m_functions.empty())
  {
    for (CFuint d = 0; d < dim; ++d)
    {
      m_functions.push_back("0.");
    }
  }

  // setup the functions and variables in the VectorialFunction object
  m_vFunction.setFunctions(m_functions);
  m_vars.resize(dim);
  m_vars[XX] = "x";
  m_vars[YY] = "y";
  if (dim == DIM_3D) { m_vars[ZZ] = "z"; }
  m_vFunction.setVariables(m_vars);

  try // parse them
  {
    m_vFunction.parse();
  }
  catch (Common::ParserException& e)
  {
    CFout << e.what() << "\n";
    throw; // retrow the exception to signal the error to the user
  }
}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallAdiabaticNSPvt::setup()
{
  FVMCC_BC::setup();

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  m_Tidx = dim+1;

  m_WallVel.resize(dim);
  m_bCoord.resize(dim);

  setupWallVelocity();
}

//////////////////////////////////////////////////////////////////////////////

void NoSlipWallAdiabaticNSPvt::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);

  // compute the coordinate at the boundary interface
  m_bCoord = (innerState->getCoordinates() + ghostState->getCoordinates());
  m_bCoord *= 0.5;

  // evaluate the velocity function in that point
  m_vFunction.evaluate(m_bCoord, m_WallVel);

  // apply the conditions to the ghost node
  (*ghostState)[0] = (*innerState)[0];
  for (CFuint d = 1; d < m_Tidx; ++d)
  {
    (*ghostState)[d] = 2*m_WallVel[d-1] - (*innerState)[d];
  }
  (*ghostState)[m_Tidx] = (*innerState)[m_Tidx];
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
