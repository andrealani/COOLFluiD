#include "NavierStokes/NSTurbTerm.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

template <class CVARSET, class DVARSET>
void NoSlipWallIsothermalNSGReKOPvt<CVARSET,DVARSET>::defineConfigOptions(Config::OptionList& options)
{
   options.template addConfigOption< CFreal >("yWallVelocity","Y-component of a velocity vector of the wall.");
   options.template addConfigOption< CFreal >("xWallVelocity","X-component of a velocity vector of the wall.");
   options.template addConfigOption< CFreal >("zWallVelocity","Z-component of a velocity vector of the wall.");
   options.template addConfigOption< CFreal >("CoeffMove","coefficient for the ghost node movement");
   options.template addConfigOption< CFreal >("KWall","Wall value for turbulent intensity");
   options.template addConfigOption< CFreal >("TWall","Wall value for the temperature");
   options.template addConfigOption< CFreal >("KGhostMin","Minumum turb. intensity in the ghost state");
}

//////////////////////////////////////////////////////////////////////////////

template <class CVARSET, class DVARSET>
NoSlipWallIsothermalNSGReKOPvt<CVARSET, DVARSET>::NoSlipWallIsothermalNSGReKOPvt(const std::string& name) :
  FVMCC_BC(name),
  _varSetTurb(CFNULL),
  _diffVarTurb(CFNULL),
  _xWallVelocity(0.),
  _yWallVelocity(0.),
  _zWallVelocity(0.),
  m_t(0.0),
  m_drXiXw(0.0),
  m_drXiXg(0.0),
  m_innerNode(CFNULL),
  m_tempNode(),
  m_midNode(),
  m_tempGhostNode(),
  m_faceNormal(),
  m_factor(0.0),
  m_kID(0),
  m_ghostK(0.),
  m_ghostGa(0.),
  m_ghostRe(0.)
{
  addConfigOptionsTo(this);

  setParameter("xWallVelocity",&_xWallVelocity);
  setParameter("yWallVelocity",&_yWallVelocity);
  setParameter("zWallVelocity",&_zWallVelocity);
  
  m_wallK = 0.;
  setParameter("KWall",&m_wallK);

  m_wallT = 0.;
  setParameter("TWall",&m_wallT);
  
  m_ghostKMin = 0.0;
  setParameter("KGhostMin",&m_ghostKMin);
  
  m_coeff = 2.0;
  setParameter("CoeffMove",&m_coeff);
}

//////////////////////////////////////////////////////////////////////////////

template <class CVARSET, class DVARSET>
NoSlipWallIsothermalNSGReKOPvt<CVARSET, DVARSET>::~NoSlipWallIsothermalNSGReKOPvt()
{
}

//////////////////////////////////////////////////////////////////////////////

template <class CVARSET, class DVARSET>
void NoSlipWallIsothermalNSGReKOPvt<CVARSET, DVARSET>::setup()
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::MathTools;

  CFAUTOTRACE;

  FVMCC_BC::setup();

  _varSetTurb = getMethodData().getUpdateVar().template d_castTo<CVARSET>();
  _diffVarTurb = getMethodData().getDiffusiveVar().template d_castTo<DVARSET>();

  _xWallVelocity /= _varSetTurb->getModel()->getVelRef();
  _yWallVelocity /= _varSetTurb->getModel()->getVelRef();
  _zWallVelocity /= _varSetTurb->getModel()->getVelRef();

  m_wallT /= _varSetTurb->getModel()->getTempRef();

  const CFuint dim = Framework::PhysicalModelStack::getActive()->getDim();
  m_tempNode.resize(dim);
  m_midNode.resize(dim);
  m_tempGhostNode.resize(dim);
  m_faceNormal.resize(dim);

  cf_assert(m_wallK >= 0.0);

  // the ID of the k variable
  m_kID = dim + 2;
}

//////////////////////////////////////////////////////////////////////////////

template <class CVARSET, class DVARSET>
void NoSlipWallIsothermalNSGReKOPvt<CVARSET, DVARSET>::setGhostState(Framework::GeometricEntity *const face)
{
 using namespace std;
 using namespace COOLFluiD::Framework;
 using namespace COOLFluiD::Common;
 using namespace COOLFluiD::MathTools;

  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const CFuint faceID = face->getID();
  const CFuint startID = faceID*dim;

  DataHandle<CFreal> normals = socket_normals.getDataHandle();

  // set the current normal
  for (CFuint i = 0; i < dim; ++i) {
    m_faceNormal[i] = normals[startID + i];
  }

  // compute the original position of the ghost state @see ComputeDummyState
  const Node& firstNode = *face->getNode(0);
  const CFreal k = - MathFunctions::innerProd(m_faceNormal, firstNode);
  const CFreal n2 = MathFunctions::innerProd(m_faceNormal, m_faceNormal);
  cf_assert(std::abs(n2) > 0.0);
  State *const innerState = face->getState(0);
  m_innerNode = &innerState->getCoordinates();
  m_t = (MathFunctions::innerProd(m_faceNormal,*m_innerNode) + k)/n2;
  m_tempGhostNode = (*m_innerNode) - 2.*m_t*m_faceNormal;

  // this middle node is by construction on the boundary face
  m_midNode = 0.5*(*m_innerNode + m_tempGhostNode);

  // first calculate the "unmodified distances" inner-wall, inner-ghost
  m_drXiXw = MathTools::MathFunctions::getDistance(*m_innerNode,m_midNode);
  m_drXiXg = MathTools::MathFunctions::getDistance(*m_innerNode, m_tempGhostNode);

  setGhostStateImpl(*innerState, *face->getState(1));
}


//////////////////////////////////////////////////////////////////////////////

template <class CVARSET, class DVARSET>
void NoSlipWallIsothermalNSGReKOPvt<CVARSET, DVARSET>::setGhostStateImpl
(const Framework::State& innerState, Framework::State& ghostState)
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::MathTools;

  // here a fix is needed in order to have always m_ghostK > 0
  // dynamic relocation of the ghost state: the position of the
  // ghost state is locally changed, and the BC is imposed
  // using a weighted average of ghost state (in the new location)
  // and inner state
  const CFreal innerK = innerState[m_kID];
  repositionNode(innerK, m_ghostK);
  
  const CFreal innerP = innerState[0];
  const CFuint TID = m_kID - 1;
  const CFuint omegaID = m_kID + 1;
  const CFuint GaID    = m_kID + 2;
  const CFuint ReID    = m_kID + 3;
  const CFreal innerT  = innerState[TID];
  const CFreal innerDensity = innerP / (_varSetTurb->getModel()->getR() * innerT);
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  CFreal wallPressure = innerP;
  cf_assert(innerP > 0.);
  cf_assert(innerT > 0.);

  // reset the ghost node with the new position
  ghostState.getCoordinates() = m_tempGhostNode;

  // here we assume that we are in PUVTGReKOega
  ghostState[0] = wallPressure;
  linearInterpolate(innerState[1], _xWallVelocity, ghostState[1]);
  linearInterpolate(innerState[2], _yWallVelocity, ghostState[2]);
  if (dim == DIM_3D) {
    linearInterpolate(innerState[3], _zWallVelocity, ghostState[3]);
  }
  
  // interpolate the wall temperature
  linearInterpolate(innerState[TID], m_wallT, ghostState[TID]);
  ghostState[m_kID] = m_ghostK;
  
  const CFuint nbTurbVars = _varSetTurb->getModel()->getNbScalarVars(0);
   //cout << "NOSPLI" << nbTurbVars << endl;
     if(nbTurbVars >= 2){

    //Compute distance to innerstate
      CFreal y0 = m_drXiXw;

      //avoid too small distances
      y0 = max(y0, 10.e-10);

      const CFreal pdim =  innerState[0] * _varSetTurb->getModel()->getPressRef();
      const CFreal Tdim =  innerState[TID] * _varSetTurb->getModel()->getTempRef();

      const CFreal mu = _diffVarTurb->getModel().getDynViscosityDim(pdim, Tdim)/
	(_diffVarTurb->getModel().getReferencePhysicalData())[Physics::NavierStokes::NSTurbTerm::MU];

      CFreal nu = mu / innerDensity;

      //this is not the best, but it avoids having to code another BC! because I
      //would have to dynamic cast to the GReKOega varset to get the beta1
      const CFreal beta1 = 0.075;

      ///@todo here should this be adimensionalized (by the distance)???
      //Menter's definition
      const CFreal omegaWall = (10. * 6. * nu) / (beta1 * y0 * y0);

      //Wilcox's definition
      //CFreal omegaWall = (1. * 6. * nu) / (beta1 * y0 * y0);
      linearInterpolate(innerState[omegaID], omegaWall, ghostState[omegaID]);

      if(ghostState[omegaID] < 0.) {
	ghostState[omegaID] = innerState[omegaID];
      }

           
     const CFreal GaWall = 1e-7;
     //BC for Ga: zero flux at the Wall 
    //const CFreal innerGa = innerState[GaID];      
     //repositionNode(innerGa, m_ghostGa); 
    // ghostState[GaID]= m_ghostGa;
     //ghostState[GaID]= innerState[GaID];
      linearInterpolate(innerState[GaID], GaWall, ghostState[GaID]);
   // cout << "Ga" << ghostState[GaID]<< endl;
      
     const CFreal ReWall = 1e-7;
     //BC for Re: zero flux at the Wall 
    //const CFreal innerRe = innerState[GaID];      
     //repositionNode(innerRe, m_ghostRe); 
     //ghostState[ReID]= m_ghostRe;
    // ghostState[ReID]= innerState[ReID];
     linearInterpolate(innerState[ReID], ReWall, ghostState[ReID]);
   // cout << "Re" << ghostState[ReID]<< endl;

  }
}

//////////////////////////////////////////////////////////////////////////////

template <class CVARSET, class DVARSET>
void NoSlipWallIsothermalNSGReKOPvt<CVARSET,DVARSET>::repositionNode
(const CFreal& innerValue, CFreal& ghostValue)
{
  using namespace std;
  using namespace COOLFluiD::Framework;
  using namespace COOLFluiD::Common;
  using namespace COOLFluiD::MathTools;

  ghostValue = 2.*m_wallK - innerValue;

  m_factor = 1.0;

  while (!(ghostValue >= m_ghostKMin)) {
    // cout << "in " << ghostValue << endl;
    m_factor *= m_coeff;

    // new position of the ghost node
    m_tempNode = ((m_factor - 1.)*m_midNode  + m_tempGhostNode)/m_factor;

    m_drXiXg = MathTools::MathFunctions::getDistance(*m_innerNode, m_tempNode);

    // new temperature in the ghost state
    linearInterpolate(innerValue, m_wallK, ghostValue);

    // cout << "innerValue  = " << innerValue << endl;
    // cout << "ghostValue  = " << ghostValue << endl;
    // cout << "m_drXiXg = " << m_drXiXg << endl;
    // cout << "m_drXiXw = " << m_drXiXw << endl;

    // move the ghost to the new position
    m_tempGhostNode = m_tempNode;
  }

  // cout << "out ghost = " << ghostValue << endl << endl;
  // cout << "out inner = " << innerValue << endl << endl;

  if (ghostValue < 0.0) cout << "ghostValue < 0 => " << ghostValue << endl;
  cf_assert(ghostValue >= 0.);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
