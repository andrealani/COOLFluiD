//#include "FluxReconstructionMethod/TetraFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionTetraP4.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

CFuint FluxReconstructionBaseFunctionTetraP4::_interpolatorID = 0;

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionBaseFunctionTetraP4::FluxReconstructionBaseFunctionTetraP4()
{
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBaseFunctionTetraP4::computeFaceJacobianDeterminant(
        const std::vector<RealVector>& mappedCoord,
        const std::vector<Framework::Node*>& nodes,
        const Framework::IntegratorPattern& pattern,
              std::vector<RealVector>& faceJacobian)
{
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionTetraP4::computeMappedCoordinates(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionTetraP4::computeMappedCoordinatesPlus1D(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
