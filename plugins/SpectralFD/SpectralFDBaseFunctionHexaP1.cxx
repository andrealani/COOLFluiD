#include "SpectralFD/HexaSpectralFDElementData.hh"
#include "SpectralFD/SpectralFDBaseFunctionHexaP1.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

CFuint SpectralFDBaseFunctionHexaP1::_interpolatorID = 0;
RealVector SpectralFDBaseFunctionHexaP1::m_ksiFac = RealVector(2);
RealVector SpectralFDBaseFunctionHexaP1::m_etaFac = RealVector(2);
RealVector SpectralFDBaseFunctionHexaP1::m_ztaFac = RealVector(2);
RealVector SpectralFDBaseFunctionHexaP1::m_solPnts1D = RealVector(2);

//////////////////////////////////////////////////////////////////////////////

SpectralFDBaseFunctionHexaP1::SpectralFDBaseFunctionHexaP1()
{
  SpectralFDElementData* sdElemData = new HexaSpectralFDElementData(getInterpolatorOrder());

  Common::SafePtr< std::vector< CFreal > > solPnts1D = sdElemData->getSolPntsLocalCoord1D();

  const CFuint nbrSolPnts = solPnts1D->size();
  cf_assert(nbrSolPnts == m_solPnts1D.size());
  for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
  {
    m_solPnts1D[iSol] = (*solPnts1D)[iSol];
  }

  delete sdElemData;
}

//////////////////////////////////////////////////////////////////////////////

void SpectralFDBaseFunctionHexaP1::computeFaceJacobianDeterminant(
        const std::vector<RealVector>& mappedCoord,
        const std::vector<Framework::Node*>& nodes,
        const Framework::IntegratorPattern& pattern,
              std::vector<RealVector>& faceJacobian)
{
    throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite difference base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector SpectralFDBaseFunctionHexaP1::computeMappedCoordinates(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite difference base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector SpectralFDBaseFunctionHexaP1::computeMappedCoordinatesPlus1D(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite difference base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
