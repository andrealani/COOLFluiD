#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionHexaP0.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

CFuint FluxReconstructionBaseFunctionHexaP0::_interpolatorID = 0;
RealVector FluxReconstructionBaseFunctionHexaP0::m_ksiFac = RealVector(1);
RealVector FluxReconstructionBaseFunctionHexaP0::m_etaFac = RealVector(1);
RealVector FluxReconstructionBaseFunctionHexaP0::m_ztaFac = RealVector(1);
RealVector FluxReconstructionBaseFunctionHexaP0::m_solPnts1D = RealVector(1);

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionBaseFunctionHexaP0::FluxReconstructionBaseFunctionHexaP0()
{
  FluxReconstructionElementData* frElemData = new HexaFluxReconstructionElementData(getInterpolatorOrder());

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

void FluxReconstructionBaseFunctionHexaP0::computeFaceJacobianDeterminant(
        const std::vector<RealVector>& mappedCoord,
        const std::vector<Framework::Node*>& nodes,
        const Framework::IntegratorPattern& pattern,
              std::vector<RealVector>& faceJacobian)
{
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionHexaP0::computeMappedCoordinates(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionHexaP0::computeMappedCoordinatesPlus1D(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
