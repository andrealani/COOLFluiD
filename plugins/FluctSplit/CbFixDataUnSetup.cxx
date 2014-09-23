#include "Framework/MethodCommandProvider.hh"

#include "FluctSplit/FluctSplit.hh"
#include "FluctSplit/CbFixDataUnSetup.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<CbFixDataUnSetup,
                      FluctuationSplitData,
                      FluctSplitModule>
aCbFixDataUnSetupProvider("CbFixDataUnSetup");

//////////////////////////////////////////////////////////////////////////////

CbFixDataUnSetup::CbFixDataUnSetup(const std::string& name) :
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
}

//////////////////////////////////////////////////////////////////////////////

CbFixDataUnSetup::~CbFixDataUnSetup()
{
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
CbFixDataUnSetup::needsSockets()
{
  std::cout << "CbFixDataUnSetup::needsSockets()\n";
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  
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

void CbFixDataUnSetup::execute()
{
  CFAUTOTRACE;

//   std::cout << "void CbFixDataUnSetup::execute()\n";

  DataHandle<CFreal> ArtViscCoeff = socket_ArtViscCoeff.getDataHandle();
  ArtViscCoeff.resize(0);
  
  DataHandle<CFreal> ViscCoeff = socket_ViscCoeff.getDataHandle();
  ViscCoeff.resize(0);

  DataHandle<CFreal> fix_active = socket_fix_active.getDataHandle();
  fix_active.resize(0);

  DataHandle<CFreal> uCsi = socket_uCsi.getDataHandle();
  uCsi.resize(0);

  DataHandle<CFreal> uEta = socket_uEta.getDataHandle();
  uEta.resize(0);

  DataHandle<CFreal> duCsidCsi = socket_duCsidCsi.getDataHandle();
  duCsidCsi.resize(0);

  DataHandle<CFreal> duEtadCsi = socket_duEtadCsi.getDataHandle();
  duEtadCsi.resize(0);

  DataHandle<CFreal> duCsidEta = socket_duCsidEta.getDataHandle();
  duCsidEta.resize(0);

  DataHandle<CFreal> duEtadEta = socket_duEtadEta.getDataHandle();
  duEtadEta.resize(0);
  
  DataHandle<CFreal> dpdCsi = socket_dpdCsi.getDataHandle();
  dpdCsi.resize(0);
  
  DataHandle<CFreal> dpdEta = socket_dpdEta.getDataHandle();
  dpdEta.resize(0);  
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
