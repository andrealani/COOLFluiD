#include "SpectralFV/SpectralFVBaseFunctionFaceTriagP1.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFV {

//////////////////////////////////////////////////////////////////////////////

CFuint SpectralFVBaseFunctionFaceTriagP1::_interpolatorID = 0;

//////////////////////////////////////////////////////////////////////////////

SpectralFVBaseFunctionFaceTriagP1::SpectralFVBaseFunctionFaceTriagP1()
{
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFVBaseFunctionFaceTriagP1::computeFaceJacobianDeterminant(
        const std::vector<RealVector>& m_mappedCoord,
        const std::vector<Framework::Node*>& nodes,
        const Framework::IntegratorPattern& pattern,
              std::vector<RealVector>& faceJacobian)
{
    throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite volume base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector SpectralFVBaseFunctionFaceTriagP1::computeMappedCoordinates(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite volume base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector SpectralFVBaseFunctionFaceTriagP1::computeMappedCoordinatesPlus1D(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite volume base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFV

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
