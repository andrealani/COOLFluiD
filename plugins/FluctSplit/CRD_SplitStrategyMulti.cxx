#include "FluctSplit/CRD_SplitStrategyMulti.hh"

#include "Framework/MethodStrategyProvider.hh"
#include "Framework/ContourIntegrator.hh"
#include "FluctSplit/FluctSplit.hh"
#include "Common/NotImplementedException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

MethodStrategyProvider<CRD_SplitStrategyMulti,
                       FluctuationSplitData,
                       FluctuationSplitStrategy,
                       FluctSplitModule>
crdFluctSplitStrategyMultiProvider("CRDMulti");

//////////////////////////////////////////////////////////////////////////////

CRD_SplitStrategyMulti::CRD_SplitStrategyMulti(const std::string& name) :
  CRD_SplitStrategy(name),
  m_systemSplitter(CFNULL),
  m_scalarSplitter(CFNULL)  
{ 
}

//////////////////////////////////////////////////////////////////////////////

CRD_SplitStrategyMulti::~CRD_SplitStrategyMulti()
{
}

//////////////////////////////////////////////////////////////////////////////

void CRD_SplitStrategyMulti::setup()
{
  cf_assert(getMethodData().isMultipleSplitter());
  
  m_scalarSplitter = getMethodData().getScalarSplitter();
  m_systemSplitter = getMethodData().getSysSplitter();
    
  // call parent method after (otherwise some common data will 
  // not be initialized correctly)
  CRD_SplitStrategy::setup();
}

//////////////////////////////////////////////////////////////////////////////

void CRD_SplitStrategyMulti::computeFluctuation(vector<RealVector>& residual)
{
 //  setCurrentCell();

//   DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
//   DistributionData& ddata = getMethodData().getDistributionData();
   
//   // compute the residual and the upwind parameters k in this cell
//   if (!getMethodData().isScalarFirst()) {
//     m_systemSplitter->computeK(*ddata.states,normals[ddata.cellID]);
//     m_scalarSplitter->computeK(*ddata.states,normals[ddata.cellID]);
//   }
//   else {
//      m_scalarSplitter->computeK(*ddata.states,normals[ddata.cellID]); 
//      m_systemSplitter->computeK(*ddata.states,normals[ddata.cellID]);
//   }
  
//   computeFluxIntegral();
  
//   // transform fluxes + source term to distribute variables
//   ddata.phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(&m_phiT);
  
//   // distribute the residual 
//   if (!getMethodData().isScalarFirst()) {
//     m_systemSplitter->distributePart(residual);
//     m_scalarSplitter->distributePart(residual);
//   }
//   else {
//     m_scalarSplitter->distributePart(residual); 
//     m_systemSplitter->distributePart(residual);
//   }

  setCurrentCell();
  
  DataHandle< InwardNormalsData*> normals = socket_normals.getDataHandle();
  DistributionData& ddata = getMethodData().getDistributionData();

  // compute the residual and the upwind parameters k in this cell
  if (!getMethodData().isScalarFirst()) {
    m_systemSplitter->computeK(*ddata.states,normals[ddata.cellID]);
    m_scalarSplitter->computeK(*ddata.states,normals[ddata.cellID]);
  }
  else {
     m_scalarSplitter->computeK(*ddata.states,normals[ddata.cellID]); 
     m_systemSplitter->computeK(*ddata.states,normals[ddata.cellID]);
  }
  
  computeFluxIntegral();

  const CFuint nbST = getMethodData().getSourceTermSplitter()->size();
  if (nbST == 1) {
   ddata.sourceTermID = 0;
    getMethodData().getSourceTermSplitter(0)->computeSourceTerm(*normals[ddata.cellID]);

    if (getMethodData().includeSourceInFlux()) {
      // in this case converctive and source fluctuations will be distributed together
      m_phiT -= ddata.phiS;
    }

    // transform fluxes + source term to distribute variables
    ddata.phi = *getMethodData().getSolutionToDistribMatTrans()->transformFromRef(&m_phiT);
    
    // distribute the residual 
    if (!getMethodData().isScalarFirst()) {
      m_systemSplitter->distributePart(residual);
      m_scalarSplitter->distributePart(residual);
    }
    else {
      m_scalarSplitter->distributePart(residual); 
      m_systemSplitter->distributePart(residual);
    }
    getMethodData().getSourceTermSplitter(0)->distribute(residual);
  }
  else {
    throw Common::NotImplementedException(FromHere(), "CRD_SplitStrategyMulti::computeFluctuation()");
  }
  
  //  if (getMethodData().hasArtificialDiff()) addArtificialDiff( residual );
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD
