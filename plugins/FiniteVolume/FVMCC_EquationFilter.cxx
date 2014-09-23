#include "Framework/VarSetTransformer.hh"
#include "Framework/MethodStrategyProvider.hh"

#include "FiniteVolume/FiniteVolume.hh"
#include "FiniteVolume/FVMCC_EquationFilter.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////


void FVMCC_EquationFilter::defineConfigOptions(Config::OptionList& options)
{
}
      
//////////////////////////////////////////////////////////////////////////////

FVMCC_EquationFilter::FVMCC_EquationFilter(const std::string& name) :
  EquationFilter<CellCenterFVMData>(name),
  socket_normals("normals"),  
  socket_faceAreas("faceAreas"),
  socket_updateCoeff("updateCoeff"),  
  socket_isOutward("isOutward"),
  socket_nstates("nstates"),
  socket_gstates("gstates")
{
}
      
//////////////////////////////////////////////////////////////////////////////

FVMCC_EquationFilter::~FVMCC_EquationFilter()
{
}

//////////////////////////////////////////////////////////////////////////////

void FVMCC_EquationFilter::setup()
{
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_EquationFilter::unsetup()
{
}
      
//////////////////////////////////////////////////////////////////////////////

void FVMCC_EquationFilter::configure ( Config::ConfigArgs& args )
{
  EquationFilter<CellCenterFVMData>::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

vector<SafePtr<BaseDataSocketSink> > FVMCC_EquationFilter::needsSockets()
{
  vector<SafePtr<BaseDataSocketSink> > result  = EquationFilter<CellCenterFVMData>::needsSockets();
  result.push_back(&socket_normals);
  result.push_back(&socket_faceAreas);
  result.push_back(&socket_updateCoeff); 
  result.push_back(&socket_isOutward);
  result.push_back(&socket_nstates);
  result.push_back(&socket_gstates);
  
  return result;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace FiniteVolume

} // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
