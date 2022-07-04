#include "FluxReconstructionBaseFunctionTriagP6.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

CFuint FluxReconstructionBaseFunctionTriagP6::_interpolatorID = 0;
RealVector FluxReconstructionBaseFunctionTriagP6::m_ksiFac = RealVector(7);
RealVector FluxReconstructionBaseFunctionTriagP6::m_etaFac = RealVector(7);
RealVector FluxReconstructionBaseFunctionTriagP6::m_solPnts1D = RealVector(7);

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionBaseFunctionTriagP6::FluxReconstructionBaseFunctionTriagP6()
{
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBaseFunctionTriagP6::computeFaceJacobianDeterminant(
        const std::vector<RealVector>& m_mappedCoord,
        const std::vector<Framework::Node*>& nodes,
        const Framework::IntegratorPattern& pattern,
              std::vector<RealVector>& faceJacobian)
{
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionTriagP6::computeMappedCoordinates(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionTriagP6::computeMappedCoordinatesPlus1D(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

