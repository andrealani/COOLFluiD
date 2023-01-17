//#include "FluxReconstructionMethod/TetraFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionTetraP1.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

CFuint FluxReconstructionBaseFunctionTetraP1::_interpolatorID = 0;
RealVector FluxReconstructionBaseFunctionTetraP1::m_ksiFac = RealVector(2);
RealVector FluxReconstructionBaseFunctionTetraP1::m_etaFac = RealVector(2);
RealVector FluxReconstructionBaseFunctionTetraP1::m_ztaFac = RealVector(2);
RealVector FluxReconstructionBaseFunctionTetraP1::m_solPnts1D = RealVector(2);

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionBaseFunctionTetraP1::FluxReconstructionBaseFunctionTetraP1()
{
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBaseFunctionTetraP1::computeFaceJacobianDeterminant(
        const std::vector<RealVector>& mappedCoord,
        const std::vector<Framework::Node*>& nodes,
        const Framework::IntegratorPattern& pattern,
              std::vector<RealVector>& faceJacobian)
{
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionTetraP1::computeMappedCoordinates(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionTetraP1::computeMappedCoordinatesPlus1D(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
