#include "SpectralFD/HexaSpectralFDElementData.hh"
#include "SpectralFD/SpectralFDBaseFunctionHexaP4.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

CFuint SpectralFDBaseFunctionHexaP4::_interpolatorID = 0;
RealVector SpectralFDBaseFunctionHexaP4::m_ksiFac = RealVector(5);
RealVector SpectralFDBaseFunctionHexaP4::m_etaFac = RealVector(5);
RealVector SpectralFDBaseFunctionHexaP4::m_ztaFac = RealVector(5);
RealVector SpectralFDBaseFunctionHexaP4::m_solPnts1D = RealVector(5);

//////////////////////////////////////////////////////////////////////////////

SpectralFDBaseFunctionHexaP4::SpectralFDBaseFunctionHexaP4()
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

void SpectralFDBaseFunctionHexaP4::computeFaceJacobianDeterminant(
        const std::vector<RealVector>& mappedCoord,
        const std::vector<Framework::Node*>& nodes,
        const Framework::IntegratorPattern& pattern,
              std::vector<RealVector>& faceJacobian)
{
    throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite difference base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector SpectralFDBaseFunctionHexaP4::computeMappedCoordinates(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite difference base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector SpectralFDBaseFunctionHexaP4::computeMappedCoordinatesPlus1D(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite difference base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
