#include "Framework/MethodCommandProvider.hh"
#include "Framework/MeshData.hh"

#include "SpectralFV/SpectralFV.hh"
#include "SpectralFV/NullCommandSpectralFV.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullCommandSpectralFV, SpectralFVMethodData, SpectralFVModule> NullCommandSpectralFVProvider("NullCommand");

//////////////////////////////////////////////////////////////////////////////

NullCommandSpectralFV::NullCommandSpectralFV(const std::string& name) :
  SpectralFVMethodCom(name)
{
}

//////////////////////////////////////////////////////////////////////////////

NullCommandSpectralFV::~NullCommandSpectralFV()
{
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD
