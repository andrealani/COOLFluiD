#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionFaceTriagP4.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

CFuint FluxReconstructionBaseFunctionFaceTriagP4::_interpolatorID = 0;

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionBaseFunctionFaceTriagP4::FluxReconstructionBaseFunctionFaceTriagP4()
{
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBaseFunctionFaceTriagP4::computeFaceJacobianDeterminant(
        const std::vector<RealVector>& mappedCoord,
        const std::vector<Framework::Node*>& nodes,
        const Framework::IntegratorPattern& pattern,
              std::vector<RealVector>& faceJacobian)
{
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionFaceTriagP4::computeMappedCoordinates(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionFaceTriagP4::computeMappedCoordinatesPlus1D(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
