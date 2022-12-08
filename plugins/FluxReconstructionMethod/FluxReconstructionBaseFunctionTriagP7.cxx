#include "FluxReconstructionBaseFunctionTriagP7.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

CFuint FluxReconstructionBaseFunctionTriagP7::_interpolatorID = 0;
RealVector FluxReconstructionBaseFunctionTriagP7::m_ksiFac = RealVector(8);
RealVector FluxReconstructionBaseFunctionTriagP7::m_etaFac = RealVector(8);
RealVector FluxReconstructionBaseFunctionTriagP7::m_solPnts1D = RealVector(8);

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionBaseFunctionTriagP7::FluxReconstructionBaseFunctionTriagP7()
{
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBaseFunctionTriagP7::computeFaceJacobianDeterminant(
        const std::vector<RealVector>& m_mappedCoord,
        const std::vector<Framework::Node*>& nodes,
        const Framework::IntegratorPattern& pattern,
              std::vector<RealVector>& faceJacobian)
{
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionTriagP7::computeMappedCoordinates(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionTriagP7::computeMappedCoordinatesPlus1D(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

