#include "FluxReconstructionBaseFunctionTriagP3.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

CFuint FluxReconstructionBaseFunctionTriagP3::_interpolatorID = 0;
RealVector FluxReconstructionBaseFunctionTriagP3::m_ksiFac = RealVector(4);
RealVector FluxReconstructionBaseFunctionTriagP3::m_etaFac = RealVector(4);
RealVector FluxReconstructionBaseFunctionTriagP3::m_solPnts1D = RealVector(4);

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionBaseFunctionTriagP3::FluxReconstructionBaseFunctionTriagP3()
{
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBaseFunctionTriagP3::computeFaceJacobianDeterminant(
        const std::vector<RealVector>& m_mappedCoord,
        const std::vector<Framework::Node*>& nodes,
        const Framework::IntegratorPattern& pattern,
              std::vector<RealVector>& faceJacobian)
{
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionTriagP3::computeMappedCoordinates(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionTriagP3::computeMappedCoordinatesPlus1D(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
