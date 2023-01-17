//#include "FluxReconstructionMethod/TetraFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionTetraP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

CFuint FluxReconstructionBaseFunctionTetraP2::_interpolatorID = 0;
RealVector FluxReconstructionBaseFunctionTetraP2::m_ksiFac = RealVector(3);
RealVector FluxReconstructionBaseFunctionTetraP2::m_etaFac = RealVector(3);
RealVector FluxReconstructionBaseFunctionTetraP2::m_ztaFac = RealVector(3);
RealVector FluxReconstructionBaseFunctionTetraP2::m_solPnts1D = RealVector(3);

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionBaseFunctionTetraP2::FluxReconstructionBaseFunctionTetraP2()
{
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBaseFunctionTetraP2::computeFaceJacobianDeterminant(
        const std::vector<RealVector>& mappedCoord,
        const std::vector<Framework::Node*>& nodes,
        const Framework::IntegratorPattern& pattern,
              std::vector<RealVector>& faceJacobian)
{
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionTetraP2::computeMappedCoordinates(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionTetraP2::computeMappedCoordinatesPlus1D(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
