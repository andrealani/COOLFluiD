#include "Framework/MethodStrategyProvider.hh"

// #define PRINT_DEBUG(x) { x }
#define PRINT_DEBUG(x) {}

// basic algorithms
#include "FluctSplit/MetaSchemes/FluctSplitMeta.hh"
#include "FluctSplit/MetaSchemes/FSMHO.hh"
#include "FluctSplit/MetaSchemes/ElemGeo.hh"
#include "FluctSplit/MetaSchemes/Physics.hh"
#include "FluctSplit/MetaSchemes/IRD.hh"
#include "FluctSplit/MetaSchemes/LRD.hh"

// schemes
#include "FluctSplit/MetaSchemes/NschemeT.hh"

// physics
#include "FluctSplit/MetaSchemes/Euler2D.hh"
#include "FluctSplit/MetaSchemes/Euler3D.hh"

// integrators
#include "FluctSplit/MetaSchemes/CIntegralTriagP1O1.hh"
#include "FluctSplit/MetaSchemes/CIntegralTriagP2O1.hh"
#include "FluctSplit/MetaSchemes/CIntegralTriagP2O2.hh"
#include "FluctSplit/MetaSchemes/CIntegralTriagP2O3.hh"
#include "FluctSplit/MetaSchemes/CIntegralTetraP1O1.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::FluctSplit;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

// TriagP1P1_Euler2D_NScheme_LRD
typedef FSMHO< ElemGeo<TriagP1,TriagP1>,
               Physics<Euler2D>,
               LRD< NschemeT > > FSMHO_0001;

MethodStrategyProvider < FSMHO_0001,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitMetaModule>
FSMHO_0001_Provider(FSMHO_0001::getSchemeName());

//////////////////////////////////////////////////////////////////////////////

// TriagP2P2_Euler2D_NScheme_LRD
typedef FSMHO< ElemGeo<TriagP2,TriagP2>,
               Physics<Euler2D>,
               LRD< NschemeT > > FSMHO_0002;

MethodStrategyProvider < FSMHO_0002,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitMetaModule>
FSMHO_0002_Provider(FSMHO_0002::getSchemeName());

//////////////////////////////////////////////////////////////////////////////

// TriagP1P1_Euler2D_NScheme_IRDContourResidual1
typedef FSMHO< ElemGeo<TriagP1,TriagP1>,
               Physics<Euler2D>,
               IRD< NschemeT, ContourResidual<CFPolyOrder::ORDER1> > > FSMHO_0003;

MethodStrategyProvider < FSMHO_0003,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitMetaModule>
FSMHO_0003_Provider(FSMHO_0003::getSchemeName());

//////////////////////////////////////////////////////////////////////////////

// TriagP2P2_Euler2D_NScheme_IRDContourResidual1
typedef FSMHO< ElemGeo<TriagP2,TriagP2>,
               Physics<Euler2D>,
               IRD< NschemeT, ContourResidual<CFPolyOrder::ORDER1> > > FSMHO_0004;

MethodStrategyProvider < FSMHO_0004,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitMetaModule>
FSMHO_0004_Provider(FSMHO_0004::getSchemeName());

//////////////////////////////////////////////////////////////////////////////

// TriagP2P2_Euler2D_NScheme_IRDContourResidual2
typedef FSMHO< ElemGeo<TriagP2,TriagP2>,
               Physics<Euler2D>,
               IRD< NschemeT, ContourResidual<CFPolyOrder::ORDER2> > > FSMHO_0005;

MethodStrategyProvider < FSMHO_0005,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitMetaModule>
FSMHO_0005_Provider(FSMHO_0005::getSchemeName());

//////////////////////////////////////////////////////////////////////////////

// TriagP2P2_Euler2D_NScheme_IRDContourResidual3
typedef FSMHO< ElemGeo<TriagP2,TriagP2>,
               Physics<Euler2D>,
               IRD< NschemeT, ContourResidual<CFPolyOrder::ORDER3> > > FSMHO_0006;

MethodStrategyProvider < FSMHO_0006,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitMetaModule>
FSMHO_0006_Provider(FSMHO_0006::getSchemeName());

//////////////////////////////////////////////////////////////////////////////

// TetraP1P1_Euler2D_NScheme_IRDContourResidual1
typedef FSMHO< ElemGeo<TetraP1,TetraP1>,
               Physics<Euler3D>,
               IRD< NschemeT, ContourResidual<CFPolyOrder::ORDER1> > > FSMHO_0007;

MethodStrategyProvider < FSMHO_0007,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitMetaModule>
FSMHO_0007_Provider(FSMHO_0007::getSchemeName());

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit

} // namespace COOLFluiD
