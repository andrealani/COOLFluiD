#include "Framework/MethodStrategyProvider.hh"
#include "Framework/MultiScalarTerm.hh"

#include "NavierStokes/EulerTerm.hh"
#include "FluctSplit/FluctSplitNavierStokes.hh"
#include "FluctSplit/FluctuationSplitData.hh"
#include "FluctSplit/BSchemeSys.hh"
#include "FluctSplit/BSchemeCSys.hh"
#include "FluctSplit/BSchemeCSysScalar.hh"
#include "FluctSplit/BxSchemeSys.hh"
#include "FluctSplit/BxSchemeSysLocal.hh"
#include "FluctSplit/BSysPSICScalar.hh"

//#include "FluctSplit/BLDAFVSchemeCSys.hh"
#include "FluctSplit/BLDAHLLSchemeCSys.hh"

#include "FluctSplit/SUPGSchemeCSys_SC.hh"


//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////
      
MethodStrategyProvider<BxSchemeSys<BSchemeSys, EulerTerm>,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitNavierStokesModule>
aBxSchemeSysProvider("SysBx");

MethodStrategyProvider<BxSchemeSys<BSchemeCSys, EulerTerm>,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitNavierStokesModule>
aBCxSchemeSysProvider("SysBCx");

MethodStrategyProvider<BxSchemeSysLocal<BSchemeSys, EulerTerm>,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitNavierStokesModule>
aBxSchemeSysLocalProvider("SysBxLocal");

MethodStrategyProvider<BxSchemeSysLocal<BSchemeCSys, EulerTerm>,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitNavierStokesModule>
aBCxSchemeSysLocalProvider("SysBCxLocal");

MethodStrategyProvider<BxSchemeSys<BSchemeSys, MultiScalarTerm<EulerTerm> >,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitNavierStokesModule>
aBxMSSchemeSysProvider("SysBxMS");
      
MethodStrategyProvider<BxSchemeSys<BSchemeCSys, MultiScalarTerm<EulerTerm> >,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitNavierStokesModule>
aBCxMSSchemeSysProvider("SysBCxMS");
    
MethodStrategyProvider<BxSchemeSysLocal<BSchemeSys, MultiScalarTerm<EulerTerm> >,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitNavierStokesModule>
aBxMSSchemeSysLocalProvider("SysBxMSLocal");

MethodStrategyProvider<BxSchemeSysLocal<BSchemeCSys, MultiScalarTerm<EulerTerm> >,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitNavierStokesModule>
aBCxMSSchemeSysLocalProvider("SysBCxMSLocalMS");

MethodStrategyProvider<BxSchemeSys<BSysPSICScalar, EulerTerm>,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitNavierStokesModule>
aBxSysPSIScalarSchemeSysProvider("BxSysPSICScalar");

MethodStrategyProvider<BxSchemeSys<BSysPSICScalar, MultiScalarTerm<EulerTerm> >,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitNavierStokesModule>
aBxSysPSIScalarMSSchemeSysProvider("BxSysPSICScalarMS");

MethodStrategyProvider<BxSchemeSys<BSchemeCSysScalar, EulerTerm>,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitNavierStokesModule>
aBCxSchemeSysScalarProvider("SysScalarBCx");

MethodStrategyProvider<BxSchemeSys<BSchemeCSysScalar, MultiScalarTerm<EulerTerm> >,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitNavierStokesModule>
aBCxMSSchemeSysScalarProvider("SysScalarBCxMS");

MethodStrategyProvider<BxSchemeSys<BLDAHLLSchemeCSys, EulerTerm>,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitNavierStokesModule>
aBCxLDAHLLSchemeSysProvider("SysBLDAHLLCx");

MethodStrategyProvider<BxSchemeSys<SUPGSchemeCSys_SC, EulerTerm>,
		       FluctuationSplitData,
		       Splitter,
                       FluctSplitNavierStokesModule>
aSUPG_SCxSchemeSysProvider("SysSUPG_SCx");



//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
