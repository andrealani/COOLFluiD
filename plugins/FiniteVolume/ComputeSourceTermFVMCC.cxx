#include "FiniteVolume/ComputeSourceTermFVMCC.hh"
#include "Framework/MeshData.hh"
#include "Framework/BaseTerm.hh"
#include "FiniteVolume/CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

void ComputeSourceTermFVMCC::defineConfigOptions
(Config::OptionList& options)
{
  options.addConfigOption< bool > 
    ("UseGradientLS", "Use the gradient computed with the least square reconstruction");
}
      
//////////////////////////////////////////////////////////////////////////////

ComputeSourceTermFVMCC::ComputeSourceTermFVMCC(const string& name) :
  ComputeSourceTerm<CellCenterFVMData>(name),
  socket_volumes("volumes"),
  socket_normals("normals"),
  socket_states("states"),
  socket_isOutward("isOutward"),  
  socket_uX("uX", false),
  socket_uY("uY", false),
  socket_uZ("uZ", false),
  m_ux(CFNULL),
  m_uy(CFNULL),
  m_uz(CFNULL),
  m_gradientsExist(false),
  m_velIDs(),
  m_pdataArray()
{
  addConfigOptionsTo(this);
  
  m_useGradientLS = false;
  setParameter("UseGradientLS", &m_useGradientLS);
}
      
//////////////////////////////////////////////////////////////////////////////

ComputeSourceTermFVMCC::~ComputeSourceTermFVMCC()
{
}

//////////////////////////////////////////////////////////////////////////////
  
void ComputeSourceTermFVMCC::configure ( Config::ConfigArgs& args )
{
  ComputeSourceTerm<CellCenterFVMData>::configure(args);
}

//////////////////////////////////////////////////////////////////////////////
  
void ComputeSourceTermFVMCC::setup()
{
  ComputeSourceTerm<CellCenterFVMData>::setup();
  
  // set the IDs corresponding to the velocity components
  getMethodData().getUpdateVar()->setStateVelocityIDs(m_velIDs);
  
  PhysicalModelStack::getActive()->getImplementor()->
    getConvectiveTerm()->resizePhysicalData(m_pdataArray);
  
  const CFuint dim = PhysicalModelStack::getActive()->getDim();
  
  const string namespaceName = MeshDataStack::getActive()->getPrimaryNamespace();
  const string uxName = namespaceName + "_uX";
  const string uyName = namespaceName + "_uY";
  const string uzName = namespaceName + "_uZ";
  
  const bool uxExists = MeshDataStack::getActive()->getDataStorage()->checkData(uxName);
  const bool uyExists = MeshDataStack::getActive()->getDataStorage()->checkData(uyName);
  const bool uzExists = MeshDataStack::getActive()->getDataStorage()->checkData(uzName);
  
  m_gradientsExist = (dim == DIM_1D && uxExists) || 
                     (dim == DIM_2D && uxExists && uyExists) || 
                     (dim == DIM_3D && uxExists && uyExists && uzExists);
  
  if (uxExists) {m_ux = MeshDataStack::getActive()->getDataStorage()->getData<CFreal>(uxName);}
  if (uyExists) {m_uy = MeshDataStack::getActive()->getDataStorage()->getData<CFreal>(uyName);}
  if (uzExists) {m_uz = MeshDataStack::getActive()->getDataStorage()->getData<CFreal>(uzName);}
}
      
//////////////////////////////////////////////////////////////////////////////

vector<Common::SafePtr<Framework::BaseDataSocketSink> >
ComputeSourceTermFVMCC::needsSockets()
{
  vector<Common::SafePtr<Framework::BaseDataSocketSink> > result =
    ComputeSourceTerm<CellCenterFVMData>::needsSockets();

  result.push_back(&socket_volumes);
  result.push_back(&socket_states);
  result.push_back(&socket_normals);
  result.push_back(&socket_isOutward);
  result.push_back(&socket_uX);
  result.push_back(&socket_uY);
  result.push_back(&socket_uZ);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

vector<Common::SafePtr<Framework::BaseDataSocketSource> >
ComputeSourceTermFVMCC::providesSockets()
{
  vector<Common::SafePtr<Framework::BaseDataSocketSource> > result =
    ComputeSourceTerm<CellCenterFVMData>::providesSockets();
  return result;
}
      
//////////////////////////////////////////////////////////////////////////////
      
    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
