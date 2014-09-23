#include "Framework/MethodCommandProvider.hh"

#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/CbFixDataSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<CbFixDataSetup,
                      FluctuationSplitData,
                      FluctSplitModule>
aCbFixDataSetupProvider("CbFixDataSetup");

//////////////////////////////////////////////////////////////////////////////

void CbFixDataSetup::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("MaxNbSubElems","Maximum number of sub-elements for high-order computations.");
}

//////////////////////////////////////////////////////////////////////////////

CbFixDataSetup::CbFixDataSetup(const std::string& name) :
  FluctuationSplitCom(name),
  socket_ArtViscCoeff("ArtViscCoeff"),
  socket_ViscCoeff("ViscCoeff"),
  socket_fix_active("fix_active"),
  socket_uCsi("uCsi"),
  socket_uEta("uEta"),
  socket_duCsidCsi("duCsidCsi"),
  socket_duEtadCsi("duEtadCsi"),
  socket_duCsidEta("duCsidEta"),
  socket_duEtadEta("duEtadEta"),
  socket_dpdCsi("dpdCsi"),
  socket_dpdEta("dpdEta")
{
  addConfigOptionsTo(this);
  m_maxsubelems = 1;
  setParameter("MaxNbSubElems",&m_maxsubelems);
}

//////////////////////////////////////////////////////////////////////////////

CbFixDataSetup::~CbFixDataSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
CbFixDataSetup::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result;
  
  result.push_back(&socket_ArtViscCoeff);
  result.push_back(&socket_ViscCoeff);
  
  result.push_back(&socket_fix_active);

  result.push_back(&socket_uCsi);
  result.push_back(&socket_uEta);

  result.push_back(&socket_duCsidCsi);
  result.push_back(&socket_duEtadCsi);
  result.push_back(&socket_duCsidEta);
  result.push_back(&socket_duEtadEta);
  
  result.push_back(&socket_dpdCsi);
  result.push_back(&socket_dpdEta);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
CbFixDataSetup::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void CbFixDataSetup::configure ( Config::ConfigArgs& args )
{
  FluctuationSplitCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

void CbFixDataSetup::execute()
{
  CFAUTOTRACE;

  const CFuint nbcells =
   MeshDataStack::getActive()->getTrs("InnerCells")->getLocalNbGeoEnts();

  DataHandle<CFreal> ArtViscCoeff = socket_ArtViscCoeff.getDataHandle();
  ArtViscCoeff.resize( nbcells );
  
  DataHandle<CFreal> ViscCoeff = socket_ViscCoeff.getDataHandle();
  ViscCoeff.resize( nbcells );
  
  DataHandle<CFreal> fix_active = socket_fix_active.getDataHandle();
  fix_active.resize( nbcells );

  DataHandle<CFreal> uCsi = socket_uCsi.getDataHandle();
  uCsi.resize( nbcells );
    
  DataHandle<CFreal> uEta = socket_uEta.getDataHandle();
  uEta.resize( nbcells );

  DataHandle<CFreal> duCsidCsi = socket_duCsidCsi.getDataHandle();
  duCsidCsi.resize( nbcells );

  DataHandle<CFreal> duEtadCsi = socket_duEtadCsi.getDataHandle();
  duEtadCsi.resize( nbcells );

  DataHandle<CFreal> duCsidEta = socket_duCsidEta.getDataHandle();
  duCsidEta.resize( nbcells );

  DataHandle<CFreal> duEtadEta = socket_duEtadEta.getDataHandle();
  duEtadEta.resize( nbcells );
  
  DataHandle<CFreal> dpdCsi = socket_dpdCsi.getDataHandle();
  dpdCsi.resize( nbcells );

  DataHandle<CFreal> dpdEta = socket_dpdEta.getDataHandle();
  dpdEta.resize( nbcells );
  
  //  ArtViscCoeff.resize( nbcells * m_maxsubelems );????
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD