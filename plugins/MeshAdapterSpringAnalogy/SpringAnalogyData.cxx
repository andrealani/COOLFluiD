#include "MeshAdapterSpringAnalogy/MeshAdapterSpringAnalogy.hh"

#include "SpringAnalogyData.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/PhysicalModelImpl.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace MeshAdapterSpringAnalogy {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<NullMethodCommand<SpringAnalogyData>, SpringAnalogyData, MeshAdapterSpringAnalogyModule> nullSpringAnalogyComProvider("Null");

//////////////////////////////////////////////////////////////////////////////

void SpringAnalogyData::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFuint >("NbSteps","NbSteps for modifying the boundary.");
}

//////////////////////////////////////////////////////////////////////////////

SpringAnalogyData::SpringAnalogyData(Common::SafePtr<Framework::Method> owner)
  : MeshAdapterData(owner),
    _geoWithNodesBuilder()
{
   addConfigOptionsTo(this);

  _totalNbSteps = 1;
   setParameter("NbSteps",&_totalNbSteps);

}

//////////////////////////////////////////////////////////////////////////////

SpringAnalogyData::~SpringAnalogyData()
{
}

//////////////////////////////////////////////////////////////////////////////

void SpringAnalogyData::configure ( Config::ConfigArgs& args )
{
  MeshAdapterData::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MeshAdapterSpringAnalogy

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

