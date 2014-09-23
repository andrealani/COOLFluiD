#include "SpectralFD/HexaSpectralFDElementData.hh"
#include "SpectralFD/SpectralFDBaseFunctionHexaP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

CFuint SpectralFDBaseFunctionHexaP2::_interpolatorID = 0;
RealVector SpectralFDBaseFunctionHexaP2::m_ksiFac = RealVector(3);
RealVector SpectralFDBaseFunctionHexaP2::m_etaFac = RealVector(3);
RealVector SpectralFDBaseFunctionHexaP2::m_ztaFac = RealVector(3);
RealVector SpectralFDBaseFunctionHexaP2::m_solPnts1D = RealVector(3);

//////////////////////////////////////////////////////////////////////////////

SpectralFDBaseFunctionHexaP2::SpectralFDBaseFunctionHexaP2()
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

void SpectralFDBaseFunctionHexaP2::computeFaceJacobianDeterminant(
        const std::vector<RealVector>& mappedCoord,
        const std::vector<Framework::Node*>& nodes,
        const Framework::IntegratorPattern& pattern,
              std::vector<RealVector>& faceJacobian)
{
    throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite difference base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector SpectralFDBaseFunctionHexaP2::computeMappedCoordinates(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite difference base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector SpectralFDBaseFunctionHexaP2::computeMappedCoordinatesPlus1D(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite difference base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
