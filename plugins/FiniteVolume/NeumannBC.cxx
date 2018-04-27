#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/NeumannBC.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/NamespaceSwitcher.hh"
#include "Config/PositiveLessThanOne.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Config;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NeumannBC, CellCenterFVMData, FiniteVolumeModule> 
neumannBCFVMCCProvider("NeumannBCFVMCC");
      
//////////////////////////////////////////////////////////////////////////////

void NeumannBC::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >
    ("OnlyRadialGradient", "Null the gradient in theta and/or phi");
}

//////////////////////////////////////////////////////////////////////////////

NeumannBC::NeumannBC(const std::string& name) :
  SuperInlet(name),
  socket_nstates("nstates"),
  _volumeCalculator(),
  _normalCalculator(),
  _n01(),
  _n12(),
  _n23(),
  _n30(),
  _n031(),
  _n132(),
  _n023(),
  _n041(),
  _n142(),
  _n024(),
  _tmpNormal(),
  _v1(),
  _v2(),
  _v3(),
  _v4(),
  _xproj(),
  _ncoord(),
  _lcoord(),
  _rcoord(),
  _tetraCoord(),
  _pyramCoord(),
  _pyramNormals1(),
  _pyramNormals2(),
  _states(),
  _nodes()
{
  addConfigOptionsTo(this);
  
  _onlyRadialGradient = false;
  setParameter("OnlyRadialGradient",&_onlyRadialGradient);
}
      
//////////////////////////////////////////////////////////////////////////////

NeumannBC::~NeumannBC()
{
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBC::setup()
{
  SuperInlet::setup();
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();

  if (dim == DIM_2D) {
    _n01.resize(dim);
    _n12.resize(dim);
    _n23.resize(dim);
    _n30.resize(dim);
    _states.resize(4);
    _nodes.resize(4);
  }
  
  if (dim == DIM_3D) {
    _n031.resize(dim);
    _n132.resize(dim);
    _n023.resize(dim);
    _n041.resize(dim);
    _n142.resize(dim);
    _n024.resize(dim);
    _tmpNormal.resize(dim);
    _v1.resize(dim);
    _v2.resize(dim);
    _v3.resize(dim);
    _v4.resize(dim);
    _xproj.resize(dim);
    _ncoord.resize(dim);
    _lcoord.resize(dim);
    _rcoord.resize(dim);
    _tetraCoord.resize(4, dim);
    _pyramCoord.resize(5, dim);
    _pyramNormals1.resize(5, dim);
    _pyramNormals2.resize(5, dim);
    _states.resize(6);
    _nodes.resize(6);
  }
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBC::unsetup()
{
  SuperInlet::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBC::setGhostState(GeometricEntity *const face)
{
  State *const innerState = face->getState(0);
  State *const ghostState = face->getState(1);
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  const bool hasIter = (_vars.size() == dim + 1);

  // coordinate of the boundary point
  _bCoord = (innerState->getCoordinates() + ghostState->getCoordinates());
  _bCoord *= 0.5;
  
  if (!_onlyRadialGradient) {
    if (!hasIter) {
      //Evaluate the function
      _vFunction.evaluate(_bCoord,*_input);
    }
    else {
      // [x,y,z] and iteration are fed to the evaluate function
      for (CFuint i = 0; i < dim; ++i) {
	_xyzIter[i] = _bCoord[i]; 
      }
      _xyzIter[dim] = SubSystemStatusStack::getActive()->getNbIter();
      
      //Evaluate the function
      _vFunction.evaluate(_xyzIter,*_input);
    }
    
    // we assume that (U_i - U_g)/dr = f(x,y,z) => U_g = U_i - f*dr
    const CFreal dr = MathTools::MathFunctions::getDistance
      (ghostState->getCoordinates(), innerState->getCoordinates());
    
    CFLog(DEBUG_MED, "NeumannBC::setGhostState() => (*_input) = " << *_input << ", dr = " << dr  << "\n");
    *ghostState = *innerState - (*_input)*dr;
  }
  else {
    if (PhysicalModelStack::getActive()->getDim() == DIM_2D) {
      // const bool isPerturb = this->getMethodData().isPerturb();
      // if (!isPerturb)
      // set the state values (pointers) corresponding to the vertices of the control volume
      computeControlVolume2D(_states, face);
      computeGhostWithRadialGradient2D(face);
    }
    if (PhysicalModelStack::getActive()->getDim() == DIM_3D) {
      // const bool isPerturb = this->getMethodData().isPerturb();
      // if (!isPerturb)
      // set the state values (pointers) corresponding to the vertices of the control volume
      if (face->nbNodes() == 3) {
	compute2TetraVolume(_states, face);
	computeGhostWithRadialGradient3DTetra(face);
      }
      else if (face->nbNodes() == 4) {
	compute2PyramVolume(_states, face);
	computeGhostWithRadialGradient3DPyram(face);
      }
    }
  }
    
  CFLog(DEBUG_MED, "NeumannBC::setGhostState() => ghostState = " << *ghostState << "\n");
}
      
//////////////////////////////////////////////////////////////////////////////
      
void NeumannBC::computeControlVolume2D
(vector<RealVector*>& states, GeometricEntity *const geo)
{
  DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  states[0] = &nodalStates[geo->getNode(0)->getLocalID()];
  states[1] = geo->getState(0);
  states[2] = &nodalStates[geo->getNode(1)->getLocalID()];
  states[3] = geo->getState(1);

  // construct the quadrilateral volume and the normals around the face
  Node& node0 = *geo->getNode(0);
  Node& node1 = geo->getState(0)->getCoordinates();
  Node& node2 = *geo->getNode(1);
  Node& node3 = geo->getState(1)->getCoordinates();

  _n01[XX] = -node0[YY] + node1[YY];
  _n01[YY] = -node1[XX] + node0[XX];

  _n12[XX] = -node1[YY] + node2[YY];
  _n12[YY] = -node2[XX] + node1[XX];

  _n23[XX] = -node2[YY] + node3[YY];
  _n23[YY] = -node3[XX] + node2[XX];

  _n30[XX] = -node3[YY] + node0[YY];
  _n30[YY] = -node0[XX] + node3[XX];

  const CFreal det = (node2[XX] - node0[XX])*(node3[YY] - node1[YY]) -
    (node3[XX] - node1[XX])*(node2[YY] - node0[YY]);
  
  if (det < 0.0) {
    _n01 *= -1.;
    _n12 *= -1.;
    _n23 *= -1.;
    _n30 *= -1.;
  }
  
  // store the nodes of the control volume
  _nodes[0] = &node0;
  _nodes[1] = &node1;
  _nodes[2] = &node2;
  _nodes[3] = &node3;
}
      
//////////////////////////////////////////////////////////////////////////////
      
vector<Common::SafePtr<Framework::BaseDataSocketSink> > NeumannBC::needsSockets()
{
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > result =
    SuperInlet::needsSockets();
  result.push_back(&socket_nstates);
  return result;
}
      
//////////////////////////////////////////////////////////////////////////////

void NeumannBC::computeGhostWithRadialGradient2D(GeometricEntity *const face)
{
  // Green Gauss is applied in the diamond volume to compute the gradients
  const CFuint nbVars = _states[0]->size();
  const CFreal x = _bCoord[XX];
  const CFreal y = _bCoord[YY];
  State *const ghostState = face->getState(1);
  
  for (CFuint i = 0; i < nbVars; ++i) {
    const CFreal v01 = (*_states[0])[i] + (*_states[1])[i];
    const CFreal v12 = (*_states[1])[i] + (*_states[2])[i];
    
    const CFreal aX = _n01[XX]*v01 + _n12[XX]*v12 +
      _n23[XX]*(*_states[2])[i] + _n30[XX]*(*_states[0])[i];
    
    const CFreal aY = _n01[YY]*v01 + _n12[YY]*v12 +
      _n23[YY]*(*_states[2])[i] + _n30[YY]*(*_states[0])[i];
    
    const CFreal den = y*(_n23[XX]+_n30[XX])-x*(_n23[YY]+_n30[YY]);
    cf_assert(std::abs(den) > 0.);
    
    // state 3 is ghost state
    (*ghostState)[i] = (-y*aX + x*aY)/den;
    

    /*CFLog(INFO, "ddTheta 0= " << -den*(*ghostState)[i]-y*aX + x*aY << "\n");
    
      const CFreal v23 = (*_states[2])[i] + (*ghostState)[i];
      const CFreal v30 = (*_states[0])[i] + (*ghostState)[i];
      const CFreal ddTheta =
      -y*(_n01[XX]*v01 + _n12[XX]*v12 + _n23[XX]*v23 + _n30[XX]*v30)
      +x*(_n01[YY]*v01 + _n12[YY]*v12 + _n23[YY]*v23 + _n30[YY]*v30);
    
    CFLog(INFO, "ddTheta 1= " << ddTheta << "\n");
    
    const CFreal ddTheta2 =
      -y*(aX + _n23[XX]*(*ghostState)[i] + _n30[XX]*(*ghostState)[i])
      +x*(aY + _n23[YY]*(*ghostState)[i] + _n30[YY]*(*ghostState)[i]);
    
    CFLog(INFO, "ddTheta 2= " << ddTheta2 << "\n");

    const CFreal ddTheta3 =
      -y*(aX + (_n23[XX] + _n30[XX])*(*ghostState)[i])
      +x*(aY + (_n23[YY] + _n30[YY])*(*ghostState)[i]);
    
      CFLog(INFO, "ddTheta 3= " << ddTheta3 << "\n");*/
  }
  
  CFLog(DEBUG_MAX, "NeumannBC::computeGhostWithRadialGradient2D() => ghostState = "
	<< *ghostState << "\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void NeumannBC::compute2TetraVolume
(vector<RealVector*>& states, GeometricEntity *const geo)
{
  DataHandle<RealVector> nodalStates = socket_nstates.getDataHandle();
  states[0] = &nodalStates[geo->getNode(0)->getLocalID()];
  states[1] = &nodalStates[geo->getNode(1)->getLocalID()];
  states[2] = &nodalStates[geo->getNode(2)->getLocalID()];
  states[3] = geo->getState(0);
  states[4] = geo->getState(1);
  
  // construct the 2 tetras volume and the normals around the face
  Node& node0 = *geo->getNode(0);
  Node& node1 = *geo->getNode(1);
  Node& node2 = *geo->getNode(2);
  Node& node3 = geo->getState(0)->getCoordinates();
  Node& node4 = geo->getState(1)->getCoordinates();
  
  // set the coordinates of the tetra 0132 in a matrix
  _tetraCoord.setRow(node0,0);
  _tetraCoord.setRow(node1,1);
  _tetraCoord.setRow(node2,2);
  _tetraCoord.setRow(node3,3);
  
  const CFreal volume0123 = _volumeCalculator.calculateTetraVolume(_tetraCoord);
  const bool changeNormalSign = (volume0123 < 0.0) ? true : false;
  CFreal factor = (!changeNormalSign) ? -0.5 : 0.5;
  
  _v1 = node3 - node0;
  _v2 = node1 - node0;
  MathFunctions::crossProd(_v1, _v2, _v3);
  _n031 = factor*_v3;
  
  _v1 = node3 - node1;
  _v2 = node2 - node1;
  MathFunctions::crossProd(_v1, _v2, _v3);
  _n132 = factor*_v3;
  
  _v1 = node2 - node0;
  _v2 = node3 - node0;
  MathFunctions::crossProd(_v1, _v2, _v3);
  _n023 = factor*_v3;
  
  // set the coordinates of the tetra 0124 in a matrix
  _tetraCoord.setRow(node4,3);
  
  // the factor is surely the opposite of the previous one
  factor *= -1;
  
  _v1 = node4 - node0;
  _v2 = node1 - node0;
  MathFunctions::crossProd(_v1, _v2, _v3);
  _n041 = factor*_v3;
  
  _v1 = node4 - node1;
  _v2 = node2 - node1;
  MathFunctions::crossProd(_v1, _v2, _v3);
  _n142 = factor*_v3;
  
  _v1 = node2 - node0;
  _v2 = node4 - node0;
  MathFunctions::crossProd(_v1, _v2, _v3);
  _n024 = factor*_v3;    
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBC::compute2PyramVolume
(vector<RealVector*>& states, GeometricEntity *const geo)
{
  DataHandle<RealVector> nodalStates = this->socket_nstates.getDataHandle();
  states[0] = &nodalStates[geo->getNode(0)->getLocalID()];
  states[1] = &nodalStates[geo->getNode(1)->getLocalID()];
  states[2] = &nodalStates[geo->getNode(2)->getLocalID()];
  states[3] = &nodalStates[geo->getNode(3)->getLocalID()];
  states[4] = geo->getState(0);
  states[5] = geo->getState(1);
  
  // construct the 2 tetras volume and the normals around the face
  Node& node0 = *geo->getNode(0);
  Node& node1 = *geo->getNode(1);
  Node& node2 = *geo->getNode(2);
  Node& node3 = *geo->getNode(3);
  Node& node4 = geo->getState(0)->getCoordinates();
  Node& node5 = geo->getState(1)->getCoordinates();
  
  // set the coordinates of the pyramid 01234 in a matrix
  _pyramCoord.setRow(node0,0);
  _pyramCoord.setRow(node1,1);
  _pyramCoord.setRow(node2,2);
  _pyramCoord.setRow(node3,3);
  _pyramCoord.setRow(node4,4);
  
  const CFreal volume01234 = _volumeCalculator.calculatePyramVolume(_pyramCoord);
  const bool changeNormalSign = (volume01234 < 0.0) ? true : false;
  CFreal factor = (!changeNormalSign) ? 1.0 : -1.0;
  
  // compute the face normals (here avoid to compute
  // the normal on the quad face !!!)
  _normalCalculator.computePyramNormals(_pyramCoord,_pyramNormals1);
  _pyramNormals1 *= factor;
  
  // set the coordinates of the pyramid 01235 in a matrix
  _pyramCoord.setRow(node5,4);
  
  // compute the face normals (here avoid to compute
  // the normal on the quad face !!!)
  _normalCalculator.computePyramNormals(_pyramCoord,_pyramNormals2);
  
  // the factor is surely the opposite of the previous one
  factor *= -1.;
  _pyramNormals2 *= factor;
}

//////////////////////////////////////////////////////////////////////////////

void NeumannBC::computeGhostWithRadialGradient3DTetra(GeometricEntity *const face)
{
  // Green Gauss is applied in the diamond volume to compute the gradients
  const CFuint nbVars = _states[0]->size();
  const CFreal x = _bCoord[XX];
  const CFreal y = _bCoord[YY];
  const CFreal z = _bCoord[ZZ];
  State *const ghostState = face->getState(1); // state 4
  
  for (CFuint i = 0; i < nbVars; ++i) {
    const CFreal s0 = (*_states[0])[i];
    const CFreal s1 = (*_states[1])[i];
    const CFreal s2 = (*_states[2])[i];
    const CFreal s3 = (*_states[3])[i];
    const CFreal s4 = (*_states[4])[i];
    
    const CFreal v031 = s0 + s3 + s1;
    const CFreal v132 = s1 + s3 + s2;
    const CFreal v023 = s0 + s2 + s3;
    
    const CFreal aX = _n031[XX]*v031 + _n132[XX]*v132 + _n023[XX]*v023 +
      _n041[XX]*(s0+s1) + _n142[XX]*(s1+s2) + _n024[XX]*(s0+s2);
    
    const CFreal aZ = _n031[ZZ]*v031 + _n132[ZZ]*v132 + _n023[ZZ]*v023 +
      _n041[ZZ]*(s0+s1) + _n142[ZZ]*(s1+s2) + _n024[ZZ]*(s0+s2);
    
    const CFreal den = z*(_n041[XX]+_n142[XX]+_n024[XX])-x*(_n041[ZZ]+_n142[ZZ]+_n024[ZZ]);
    cf_assert(std::abs(den) > 0.);
    
    // state 4 is ghost state
    (*ghostState)[i] = (-z*aX + x*aZ)/den;
  }
  
  CFLog(DEBUG_MAX, "NeumannBC::computeGhostWithRadialGradient3DTetra() => ghostState = "
	<< *ghostState << "\n");
}
      
//////////////////////////////////////////////////////////////////////////////

void NeumannBC::computeGhostWithRadialGradient3DPyram(GeometricEntity *const face)
{
  throw NotImplementedException
    (FromHere(), "NeumannBC::computeGhostWithRadialGradient3DPyram()");
}
      
//////////////////////////////////////////////////////////////////////////////
      
    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
