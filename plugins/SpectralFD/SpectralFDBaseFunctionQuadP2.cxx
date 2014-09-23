#include "SpectralFD/QuadSpectralFDElementData.hh"
#include "SpectralFD/SpectralFDBaseFunctionQuadP2.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

CFuint SpectralFDBaseFunctionQuadP2::_interpolatorID = 0;
RealVector SpectralFDBaseFunctionQuadP2::m_ksiFac = RealVector(3);
RealVector SpectralFDBaseFunctionQuadP2::m_etaFac = RealVector(3);
RealVector SpectralFDBaseFunctionQuadP2::m_solPnts1D = RealVector(3);

//////////////////////////////////////////////////////////////////////////////

SpectralFDBaseFunctionQuadP2::SpectralFDBaseFunctionQuadP2()
{
  SpectralFDElementData* sdElemData = new QuadSpectralFDElementData(getInterpolatorOrder());

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

void SpectralFDBaseFunctionQuadP2::computeFaceJacobianDeterminant(
        const std::vector<RealVector>& mappedCoord,
        const std::vector<Framework::Node*>& nodes,
        const Framework::IntegratorPattern& pattern,
              std::vector<RealVector>& faceJacobian)
{
    throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite difference base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector SpectralFDBaseFunctionQuadP2::computeMappedCoordinates(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite difference base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

RealVector SpectralFDBaseFunctionQuadP2::computeMappedCoordinatesPlus1D(const RealVector& coord,
                                    const std::vector<Framework::Node*>& nodes)
{
  throw Common::ShouldNotBeHereException (FromHere(),"Spectral finite difference base functions should not be used as geometrical shape functions.");
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace SpectralFD

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
