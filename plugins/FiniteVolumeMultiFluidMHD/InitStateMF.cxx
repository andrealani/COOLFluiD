#include "FiniteVolumeMultiFluidMHD/FiniteVolumeMultiFluidMHD.hh"
#include "FiniteVolumeMultiFluidMHD/InitStateMF.hh"
#include "Common/CFLog.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Common/BadValueException.hh"
#include "Framework/NamespaceSwitcher.hh"
#include <cmath>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<InitStateMF, CellCenterFVMData, FiniteVolumeMultiFluidMHDModule>
initStateMFProvider("InitStateMF");

//////////////////////////////////////////////////////////////////////////////

void InitStateMF::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

InitStateMF::InitStateMF(const std::string& name) :
  CellCenterFVMCom(name),
  socket_states("states")
{
  addConfigOptionsTo(this);
}

//////////////////////////////////////////////////////////////////////////////

InitStateMF::~InitStateMF()
{
}

//////////////////////////////////////////////////////////////////////////////

void InitStateMF::executeOnTrs()
{  
  CFLog(VERBOSE, "InitStateMF::executeOnTrs() => START\n");
  
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();
  
  for (CFuint i = 0; i < states.size(); ++i) {
    State& currState  = *states[i];
    const Node& coord = currState.getCoordinates();
    const CFreal x = coord[XX];
    const CFreal y = coord[YY];
    
    // set default to 0.
    currState = 0.;
    // here hardcode expression for each non-zero state variable (using C-style syntax ... for instance "^" doesn't exist -> use std::pow(...,...))
    currState[1] = 1.8e-3*std::exp(-((x-2e6)*(x-2e6))/(2*(40000*40000)));
    
    // currState[8] = ; adapt the following expression ...the same for the other
    //if(y<391959.798995,1.672621777e-27*10^(1.9810005627935e+01-7.5542888256189e-03*(y*1e-3)^1-1.6242962465733e-04*(y*1e-3)^2+4.0838641524838e-06*(y*1e-3)^3-4.3168415747324e-08*(y*1e-3)^4+2.5854049322152e-10*(y*1e-3)^5-9.3655454802867e-13*(y*1e-3)^6+2.0327629131765e-15*(y*1e-3)^7-2.4341611305037e-18*(y*1e-3)^8+1.2369208397681e-21*(y*1e-3)^9),if(y<994974.874372,1.672621777e-27*10^(-2.0102228762052e+03+2.8625543055849e+01*(y*1e-3)^1-1.7644301150933e-01*(y*1e-3)^2+6.2333706804043e-04*(y*1e-3)^3-1.3908154904903e-06*(y*1e-3)^4+2.0325362410282e-09*(y*1e-3)^5-1.9459993149366e-12*(y*1e-3)^6+1.1776404636220e-15*(y*1e-3)^7-4.0903048729410e-19*(y*1e-3)^8+6.2176222138913e-23*(y*1e-3)^9),if(y<1195979.8995,1.672621777e-27*10^(1.2331401919204e+08-1.0104153732651e+06*(y*1e-3)^1+3.6777219168941e+03*(y*1e-3)^2-7.8045908262529e+00*(y*1e-3)^3+1.0641749497571e-02*(y*1e-3)^4-9.6685985102399e-06*(y*1e-3)^5+5.8533346386708e-09*(y*1e-3)^6-2.2768792930054e-12*(y*1e-3)^7+5.1638972955402e-16*(y*1e-3)^8-5.2025782106265e-20*(y*1e-3)^9),if(y<1447236.1809,1.672621777e-27*10^(1.3588102106555e+08-9.2432598874381e+05*(y*1e-3)^1+2.7928203059235e+03*(y*1e-3)^2-4.9193949178004e+00*(y*1e-3)^3+5.5671011124037e-03*(y*1e-3)^4-4.1974906854201e-06*(y*1e-3)^5+2.1085982042780e-09*(y*1e-3)^6-6.8053018647910e-13*(y*1e-3)^7+1.2804217839025e-16*(y*1e-3)^8-1.0700691393828e-20*(y*1e-3)^9),if(y<1798994.97487,1.672621777e-27*10^(-2.4191887751355e+06+1.2260463931990e+04*(y*1e-3)^1-2.7413238688378e+01*(y*1e-3)^2+3.5450312486208e-02*(y*1e-3)^3-2.9175054710523e-05*(y*1e-3)^4+1.5813414390668e-08*(y*1e-3)^5-5.6287050759954e-12*(y*1e-3)^6+1.2634511841218e-15*(y*1e-3)^7-1.6126376586325e-19*(y*1e-3)^8+8.8266420494005e-24*(y*1e-3)^9),if(y<2000000.0,1.672621777e-27*10^(1.1342966769414e+09-5.4232931567978e+06*(y*1e-3)^1+1.1520909422659e+04*(y*1e-3)^2-1.4272423405720e+01*(y*1e-3)^3+1.1363034697371e-02*(y*1e-3)^4-6.0293646091789e-06*(y*1e-3)^5+2.1322103003493e-09*(y*1e-3)^6-4.8459195389539e-13*(y*1e-3)^7+6.4226316099586e-17*(y*1e-3)^8-3.7821767599565e-21*(y*1e-3)^9),1.672621777e-27*10^(16.583876609802246)))))))
    
    
  }
  
  CFLog(VERBOSE, "InitStateMF::executeOnTrs() => END\n");
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSink> >
InitStateMF::needsSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSink> > result;
  result.push_back(&socket_states);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

void InitStateMF::setup()
{
  CFAUTOTRACE;

  CellCenterFVMCom::setup();
}

//////////////////////////////////////////////////////////////////////////////

void InitStateMF::unsetup()
{
  CFAUTOTRACE;

  CellCenterFVMCom::unsetup();
}

//////////////////////////////////////////////////////////////////////////////

void InitStateMF::configure ( Config::ConfigArgs& args )
{
  CFAUTOTRACE;

  CellCenterFVMCom::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
