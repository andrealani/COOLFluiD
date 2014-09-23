#include "RoeFluxALEBDF2.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/MethodStrategyProvider.hh"
#include "Framework/GeometricEntity.hh"
#include "FiniteVolume/FiniteVolume.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<RoeFluxALEBDF2,
                       CellCenterFVMData,
                       FluxSplitter<CellCenterFVMData>,
                       FiniteVolumeModule>
roeFluxALEBDF2Provider("RoeFluxALEBDF2");
 
//////////////////////////////////////////////////////////////////////////////

void RoeFluxALEBDF2::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< bool >("UseShapeFunctions", "Use the shape functions");
}

//////////////////////////////////////////////////////////////////////////////

RoeFluxALEBDF2::RoeFluxALEBDF2(const std::string& name) :
  RoeFluxALE(name),
  socket_pastPastNodes("pastPastNodes"),
  socket_avNormals("avNormals"),
  _pastUnitNormal()
{
  addConfigOptionsTo(this);
  
  _useShapeFunctions = true;
  setParameter("UseShapeFunctions", &_useShapeFunctions);
}
      
//////////////////////////////////////////////////////////////////////////////

RoeFluxALEBDF2::~RoeFluxALEBDF2()
{
}

//////////////////////////////////////////////////////////////////////////////

void RoeFluxALEBDF2::setAbsEigenValues()
{
  // Compute the meshSpeed
  // and Modify the eigen values to account for mesh deformation
  const CFreal dt = SubSystemStatusStack::getActive()->getDT();
  const bool isFirstStep = (SubSystemStatusStack::getActive()->getNbIter() == 1) ? true : false;
  CellCenterFVMData& data = this->getMethodData(); 
  GeometricEntity *const geo = data.getCurrentFace();
  SafePtr<FVMCC_PolyRec> polyRec = data.getPolyReconstructor();
  const RealVector& shapeFunction = polyRec->getCurrentGeoShapeFunction(geo);
  
  const CFuint nbNodes = geo->nbNodes();
  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();
  
  DataHandle<Node*> pastNodes = socket_pastNodes.getDataHandle();
  DataHandle<Node*> pastPastNodes = socket_pastPastNodes.getDataHandle();
  DataHandle<Node*> futureNodes = socket_futureNodes.getDataHandle();
  DataHandle<CFreal> normals = socket_normals.getDataHandle();
  DataHandle<CFreal> avNormals = socket_avNormals.getDataHandle();
  
  if(!isFirstStep){
    const CFreal previousDt = SubSystemStatusStack::getActive()->getPreviousDT();
    const CFuint faceID = geo->getID();
    const CFuint startID = faceID*nbDim;
    const CFreal alpha = previousDt / dt;
    const CFreal xi = 1. / (1. + alpha);
    const CFreal coef1 = (1. + xi)/dt;
    const CFreal coef2 = xi / previousDt;

    for(CFuint iDim = 0; iDim < nbDim; iDim++){
      _pastUnitNormal[iDim] = normals[startID + iDim] - ((1.+xi) * avNormals[startID + iDim]);
      _pastUnitNormal[iDim] /= -xi;
      _tempUnitNormal[iDim] = avNormals[startID + iDim];
    }
    
    _pastUnitNormal.normalize();
    _tempUnitNormal.normalize();
    
    // Compute speed of the mesh at current quadrature point
    _vgn = 0.;
    for(CFuint iNode = 0; iNode < nbNodes; iNode++){
      const CFuint nodeID = (geo->getNode(iNode))->getLocalID();
      for(CFuint iDim = 0; iDim < nbDim; iDim++){
        const CFreal shapeF = (_useShapeFunctions) ? shapeFunction[iNode] : 1.;
	
	_vgn +=
	  shapeF * _tempUnitNormal[iDim] * coef1 * ( (*(futureNodes[nodeID]))[iDim] - (*(pastNodes[nodeID]))[iDim]);

        _vgn -=
	  shapeF * _pastUnitNormal[iDim] * coef2 * ( (*(pastNodes[nodeID]))[iDim] - (*(pastPastNodes[nodeID]))[iDim]);
      }
    }
  }
  // take care of the first step
  // if it is the first step, then use CrankNicholson
  else{
    _meshSpeed = 0.;
    // Compute speed of the mesh at current quadrature point
    for(CFuint iNode = 0; iNode < nbNodes; iNode++){
      const CFuint nodeID = (geo->getNode(iNode))->getLocalID();
      for(CFuint iDim = 0; iDim < nbDim; iDim++){
	const CFreal shapeF = (_useShapeFunctions) ? shapeFunction[iNode] : 1.;
        _meshSpeed[iDim] += shapeF * (*(futureNodes[nodeID]))[iDim];
        _meshSpeed[iDim] -= shapeF * (*(pastNodes[nodeID]))[iDim];
      }
    }
    
    _meshSpeed /= dt;

    //Compute vg*n 
    const RealVector& unitNormal = getMethodData().getUnitNormal();
    _vgn = _meshSpeed[0] * unitNormal[0];
    for(CFuint iDim = 1;iDim < nbDim ;iDim++){
      _vgn += _meshSpeed[iDim] * unitNormal[iDim];
    }
  }

  // Modify the eigen values (eValue = eValue - Vg.n)
  _eValues -= _vgn;

  //Set absolute value of eigenvalues
  _absEvalues = abs(_eValues);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > RoeFluxALEBDF2::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result = RoeFluxALE::needsSockets();
  
  result.push_back(&socket_pastPastNodes);
  result.push_back(&socket_avNormals);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
