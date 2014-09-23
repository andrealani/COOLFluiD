#include "HessianEE/HessianEE.hh"
#include "StdSetup.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace HessianEE {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<StdSetup, HessEEData, HessianEEModule> stdSetupProvider("StdSetup");

//////////////////////////////////////////////////////////////////////////////

StdSetup::StdSetup(std::string name) : HessEECom(name),
  socket_adapt_func("adapt_func"),
  socket_grad("grad"),
  socket_hessian("hessian"),
  socket_glob_metric("glob_metric"),
  socket_metric("metric"),
  socket_adapt_wght("adapt_wght"),
  socket_adapt_matrix("adapt_matrix"),
  socket_nodes("nodes")
{
}

//////////////////////////////////////////////////////////////////////////////

StdSetup::~StdSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
StdSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;

  result.push_back(&socket_adapt_func);
  result.push_back(&socket_grad);
  result.push_back(&socket_hessian);

  result.push_back(&socket_glob_metric);
  result.push_back(&socket_metric);

  result.push_back(&socket_adapt_wght);
  result.push_back(&socket_adapt_matrix);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
StdSetup::needsSockets()
{

  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_nodes);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

void StdSetup::execute()
{
  CFAUTOTRACE;

  DataHandle < Framework::Node*, Framework::GLOBAL > nodes = socket_nodes.getDataHandle();

  //const CFuint nbStates = states.size();
  const CFuint nbNodes = nodes.size();

  const CFuint nbDim = PhysicalModelStack::getActive()->getDim();

  RealMatrix initMtx( nbDim, nbDim, 0.0);
  RealVector initVec( 0.0, nbDim);

  DataHandle<CFreal> adapt_func = socket_adapt_func.getDataHandle();
  adapt_func.resize(nbNodes);
  adapt_func = 0.0;

  DataHandle<RealVector> grad = socket_grad.getDataHandle();
  grad.resize(nbNodes);

  DataHandle<RealMatrix> hessian = socket_hessian.getDataHandle();
  hessian.resize(nbNodes);



  DataHandle<RealMatrix> gmetric = socket_glob_metric.getDataHandle();
  gmetric.resize(nbNodes);

  DataHandle<RealMatrix> metric = socket_metric.getDataHandle();
  metric.resize(nbNodes);



  DataHandle<CFreal> adapt_wght = socket_adapt_wght.getDataHandle();
  adapt_wght.resize(nbNodes);
  adapt_wght = 0.0;

  DataHandle<RealMatrix> adapt_matrix = socket_adapt_matrix.getDataHandle();
  adapt_matrix.resize(nbNodes);



  for(CFuint i = 0; i < nbNodes; ++i)
  {

    grad[i].resize(nbDim);
    gmetric[i].resize(nbDim,nbDim);
    metric[i].resize(nbDim,nbDim);
    hessian[i].resize(nbDim,nbDim);
	adapt_matrix[i].resize(nbDim,nbDim);

    grad[i] = 0.0;
    gmetric[i] = 0.0;
    metric[i] = 0.0;
    hessian[i] = 0.0;
    adapt_matrix[i] = 0.0;
  }
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace HessianEE

  } // namespace Numerics

} // namespace COOLFluiD
