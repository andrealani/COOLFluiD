//#include "FluxReconstructionMethod/TetraFluxReconstructionElementData.hh"
#include "FluxReconstructionMethod/FluxReconstructionBaseFunctionTetraP5.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

CFuint FluxReconstructionBaseFunctionTetraP5::_interpolatorID = 0;

//////////////////////////////////////////////////////////////////////////////

FluxReconstructionBaseFunctionTetraP5::FluxReconstructionBaseFunctionTetraP5()
{
}

//////////////////////////////////////////////////////////////////////////////

void FluxReconstructionBaseFunctionTetraP5::computeFaceJacobianDeterminant(
        const std::vector<RealVector>& mappedCoord,
        const std::vector<Framework::Node*>& nodes,
        const Framework::IntegratorPattern& pattern,
              std::vector<RealVector>& faceJacobian)
{
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionTetraP5::computeMappedCoordinates(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector FluxReconstructionBaseFunctionTetraP5::computeMappedCoordinatesPlus1D(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
