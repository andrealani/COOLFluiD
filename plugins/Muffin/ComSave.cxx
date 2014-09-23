
#include "Framework/MethodCommandProvider.hh"
#include "Framework/OutputFormatter.hh"
#include "Muffin/MuffinModule.hh"
#include "Muffin/ComSave.hh"

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

namespace COOLFluiD {
  namespace Muffin {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider< ComSave,MuffinData,MuffinModule > cComSaveProvider("ComSave");

//////////////////////////////////////////////////////////////////////////////

void ComSave::execute()
{
  getMethodData().log("save...");
  const bool forceWriting = true;
  std::vector< SafePtr< OutputFormatter > > m = MethodRegistry::getInstance().getAllMethods< OutputFormatter >(getMethodData().getNamespace());
  for (CFuint i=0; i<m.size(); ++i)
    if (!m[i]->isNonRootMethod() && (m[i]->isSaveNow( forceWriting ) ) ) {
      m[i]->open();
      m[i]->write();
      m[i]->close();
    }
  getMethodData().ver("save.");
}

//////////////////////////////////////////////////////////////////////////////

  }  // namespace Muffin
}  // namespace COOLFluiD

