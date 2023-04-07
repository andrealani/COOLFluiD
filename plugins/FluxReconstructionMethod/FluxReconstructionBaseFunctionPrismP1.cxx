//#include "FluxReconstructionMethod/PrismFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionPrismP1.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

CFuint FluxReconstructionBaseFunctionPrismP1::_interpolatorID = 0;
RealVector FluxReconstructionBaseFunctionPrismP1::m_ksiFac = RealVector(2);
RealVector FluxReconstructionBaseFunctionPrismP1::m_etaFac = RealVector(2);
RealVector FluxReconstructionBaseFunctionPrismP1::m_ztaFac = RealVector(2);
RealVector FluxReconstructionBaseFunctionPrismP1::m_solPnts1D = RealVector(2);

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionBaseFunctionPrismP1::FluxReconstructionBaseFunctionPrismP1()
{
  FluxReconstructionElementData* frElemData = new PrismFluxReconstructionElementData(getInterpolatorOrder());

  Common::SafePtr< std::vector< CFreal > > solPnts1D = frElemData->getSolPntsLocalCoord1D();

  const CFuint nbrSolPnts = solPnts1D->size();
  cf_assert(nbrSolPnts == m_solPnts1D.size());
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    m_solPnts1D[iSol] = (*solPnts1D)[iSol];
  }

  delete frElemData;
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBaseFunctionPrismP1::computeFaceJacobianDeterminant(
        const std::vector<RealVector>& mappedCoord,
        const std::vector<Framework::Node*>& nodes,
        const Framework::IntegratorPattern& pattern,
              std::vector<RealVector>& faceJacobian)
{
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionPrismP1::computeMappedCoordinates(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionPrismP1::computeMappedCoordinatesPlus1D(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
