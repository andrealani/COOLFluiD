#include "ComputeSpaceResidual.hh"
#include "Framework/MeshData.hh"
#include "Framework/VolumeIntegrator.hh"
#include "Framework/LSSMatrix.hh"
#include "Framework/BlockAccumulator.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteElement {

//////////////////////////////////////////////////////////////////////////////

ComputeSpaceResidual::ComputeSpaceResidual(const std::string& name) :
  FiniteElementMethodCom(name),
  m_map_femdata(),
  socket_rhs("rhs")
{
}

//////////////////////////////////////////////////////////////////////////////

ComputeSpaceResidual::~ComputeSpaceResidual()
{
}

//////////////////////////////////////////////////////////////////////////////

void ComputeSpaceResidual::unsetup()
{
  for(CFuint i = 0; i < m_map_femdata.size(); ++i)
  {
    deletePtr(m_map_femdata[i].first);
    deletePtr(m_map_femdata[i].second);
    deletePtr(m_map_femdata[i].third);
    deletePtr(m_map_femdata[i].fourth);
  }

  // then call parent method
  FiniteElementMethodCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void ComputeSpaceResidual::setup()
{
  CFAUTOTRACE;

  // first call parent method
  FiniteElementMethodCom::setup();

  // rhs storage
  DataHandle<CFreal> rhs = socket_rhs.getDataHandle();

  SafePtr<vector<ElementTypeData> > elementType = MeshDataStack::getActive()->getElementTypeData();
  const CFuint nbElemTypes = elementType->size();
  const CFuint nbEqs = PhysicalModelStack::getActive()->getNbEq();

  // loop over types since it can happen to deal with an hybrid mesh
  for (CFuint iType = 0; iType < nbElemTypes; ++iType) {

    const CFuint nbStatesInType = (*elementType)[iType].getNbStates();

    BlockAccumulator* ptr = getMethodData().getLinearSystemSolver()[0]->
      createBlockAccumulator(nbStatesInType,nbStatesInType,nbEqs);
    RealVector* vec = new RealVector(nbStatesInType*nbEqs);
    RealMatrix* mat = new RealMatrix(nbStatesInType*nbEqs,nbStatesInType*nbEqs);

    RealVector  prototype(nbEqs);
    vector<RealVector>* residual = new vector<RealVector>(nbStatesInType,prototype);

    FElemTypeData elemTypeData(ptr,mat,vec,residual);

    m_map_femdata.insert(nbStatesInType, elemTypeData);
  }

  m_map_femdata.sortKeys();
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
ComputeSpaceResidual::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;

  result.push_back(&socket_rhs);

  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteElement

  } // namespace Numerics

} // namespace COOLFluiD
