#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionPrismP0.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

CFuint FluxReconstructionBaseFunctionPrismP0::_interpolatorID = 0;
RealVector FluxReconstructionBaseFunctionPrismP0::m_ksiFac = RealVector(1);
RealVector FluxReconstructionBaseFunctionPrismP0::m_etaFac = RealVector(1);
RealVector FluxReconstructionBaseFunctionPrismP0::m_ztaFac = RealVector(1);
RealVector FluxReconstructionBaseFunctionPrismP0::m_solPnts1D = RealVector(1);

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionBaseFunctionPrismP0::FluxReconstructionBaseFunctionPrismP0()
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

void FluxReconstructionBaseFunctionPrismP0::computeFaceJacobianDeterminant(
        const std::vector<RealVector>& mappedCoord,
        const std::vector<Framework::Node*>& nodes,
        const Framework::IntegratorPattern& pattern,
              std::vector<RealVector>& faceJacobian)
{
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionPrismP0::computeMappedCoordinates(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionPrismP0::computeMappedCoordinatesPlus1D(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
