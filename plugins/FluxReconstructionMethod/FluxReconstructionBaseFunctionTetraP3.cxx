//#include "FluxReconstructionMethod/TetraFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionTetraP3.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

CFuint FluxReconstructionBaseFunctionTetraP3::_interpolatorID = 0;
RealVector FluxReconstructionBaseFunctionTetraP3::m_ksiFac = RealVector(4);
RealVector FluxReconstructionBaseFunctionTetraP3::m_etaFac = RealVector(4);
RealVector FluxReconstructionBaseFunctionTetraP3::m_ztaFac = RealVector(4);
RealVector FluxReconstructionBaseFunctionTetraP3::m_solPnts1D = RealVector(4);

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionBaseFunctionTetraP3::FluxReconstructionBaseFunctionTetraP3()
{
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBaseFunctionTetraP3::computeFaceJacobianDeterminant(
        const std::vector<RealVector>& mappedCoord,
        const std::vector<Framework::Node*>& nodes,
        const Framework::IntegratorPattern& pattern,
              std::vector<RealVector>& faceJacobian)
{
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionTetraP3::computeMappedCoordinates(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionTetraP3::computeMappedCoordinatesPlus1D(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
