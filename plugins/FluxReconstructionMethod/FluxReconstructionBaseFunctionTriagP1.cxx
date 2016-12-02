#include "FluxReconstructionBaseFunctionTriagP1.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

CFuint FluxReconstructionBaseFunctionTriagP1::_interpolatorID = 0;
RealVector FluxReconstructionBaseFunctionTriagP1::m_ksiFac = RealVector(2);
RealVector FluxReconstructionBaseFunctionTriagP1::m_etaFac = RealVector(2);
RealVector FluxReconstructionBaseFunctionTriagP1::m_solPnts1D = RealVector(2);

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionBaseFunctionTriagP1::FluxReconstructionBaseFunctionTriagP1()
{
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBaseFunctionTriagP1::computeFaceJacobianDeterminant(
        const std::vector<RealVector>& m_mappedCoord,
        const std::vector<Framework::Node*>& nodes,
        const Framework::IntegratorPattern& pattern,
              std::vector<RealVector>& faceJacobian)
{
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionTriagP1::computeMappedCoordinates(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionTriagP1::computeMappedCoordinatesPlus1D(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
