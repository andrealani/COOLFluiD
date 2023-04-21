#include "MHD/MHD.hh"
#include "MHD/MHD3DProjectionEpsDiffPrim.hh"
#include "Environment/ObjectProvider.hh"
#include "MHD/MHDProjectionEpsTerm.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace MHD {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<MHD3DProjectionEpsDiffPrim, DiffusiveVarSet, MHDModule, 2> 
mhd3DProjectionEpsDiffPrimProvider("MHD3DProjectionEpsDiffPrim");
      
//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionEpsDiffPrim::MHD3DProjectionEpsDiffPrim(const std::string& name,
						 Common::SafePtr<Framework::PhysicalModelImpl> model) :
  MHD3DProjectionDiffPrim(name, model)
{
  string names9  = "epsP";
  string names10 = "epsM";
  addVarName(names9);
  addVarName(names10);
}

//////////////////////////////////////////////////////////////////////////////

MHD3DProjectionEpsDiffPrim::~MHD3DProjectionEpsDiffPrim()
{
}

//////////////////////////////////////////////////////////////////////////////
      
} // namespace MHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
