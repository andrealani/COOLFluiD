#ifndef COOLFluiD_FluxReconstructionMethod_FluxReconstructionBaseFunctionPrismP0_hh
#define COOLFluiD_FluxReconstructionMethod_FluxReconstructionBaseFunctionPrismP0_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/StringOps.hh"
#include "ShapeFunctions/LagrangeShapeFunction.hh"
#include "Common/NotImplementedException.hh"
#include "Common/ShouldNotBeHereException.hh"
#include "MathTools/RealMatrix.hh"
#include "FluxReconstructionMethod/PrismFluxReconstructionElementData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

 namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class provides the FR base function describing the
 * representation of the solution in a P0 Prism element.
 *
 * @author Ray Vandenhoeck
 * @author Alexander Papen
 * @author Rayan Dhib
 */
class FluxReconstructionBaseFunctionPrismP0 : public ShapeFunctions::LagrangeShapeFunction {
public:

  /**
   * Get the inheritant dimensionality of the ShapeFunction
   */
  static CFuint getDimensionality()
  {
    return DIM_3D;
  }

  /**
   * Get the number of nodal base functions in the element
   */
  static CFuint getNbNodes()
  {
    return 1;
  }

  /**
   * Get the number of faces
   */
  static CFuint getNbFaces()
  {
    return 5;
  }

  /**
   * Gets the type of CFGeoShape::Type
   */
  static CFGeoShape::Type getShape()
  {
    return CFGeoShape::PRISM;
  }

 /**
  * Gets the type of Interpolator
  */
  static CFPolyForm::Type getInterpolatorType()
  {
    return CFPolyForm::FLUXRECONSTRUCTION;
  }

  /// Gets the mapped coordinates of the DOFs
  /// @param mappedCoords where to put the coordinates
  static void getStatesMappedCoordinates(std::vector<RealVector>& mappedCoords)
  {
    throw Common::NotImplementedException (FromHere(),getName() + "::getStatesMappedCoordinates()");
  }

 /**
  * Gets the Interpolator order
  */
  static CFPolyOrder::Type getInterpolatorOrder()
  {
    return CFPolyOrder::ORDER0;
  }

  /**
   * Compute the base functions corresponding to the given
   * mapped coordinates
   */
  static void computeBaseFunctions(
          const std::vector<RealVector>& mappedCoord, std::vector<RealVector>& shapeFunc)
  {
    for (CFuint ip = 0; ip < mappedCoord.size(); ++ip) {
      computeShapeFunction(mappedCoord[ip],shapeFunc[ip]);
    }
  }

  /**
   * Compute the base functions corresponding to the given
   * mapped coordinates
   */
  static void computeShapeFunction(
        const RealVector& mappedCoord, RealVector& shapeFunc)
  {
    FluxReconstructionElementData* frElemData = new PrismFluxReconstructionElementData(getInterpolatorOrder());

    Common::SafePtr< std::vector< CFreal > > solPnts1D = frElemData->getSolPntsLocalCoord1D();

    const CFuint nbrSolPnts = solPnts1D->size();
    cf_assert(nbrSolPnts == m_solPnts1D.size());
    for (CFuint iSol = 0; iSol < nbrSolPnts; ++iSol)
    {
      m_solPnts1D[iSol] = (*solPnts1D)[iSol];
    }

    delete frElemData;
    // coordinates of output points
    const CFreal ksi = mappedCoord[KSI];
    const CFreal eta = mappedCoord[ETA];
    const CFreal zta = mappedCoord[ZTA];

    // ksi factors
    for (CFuint iSol = 0; iSol < 1; ++iSol)
    {
      const CFreal ksiSol = m_solPnts1D[iSol];
      m_ksiFac[iSol] = 1.;
      for (CFuint iFac = 0; iFac < 1; ++iFac)
      {
        if (iFac != iSol)
        {
          const CFreal ksiFac = m_solPnts1D[iFac];
          m_ksiFac[iSol] *= (ksi-ksiFac)/(ksiSol-ksiFac);
        }
      }
    }

    // eta factors
    for (CFuint iSol = 0; iSol < 1; ++iSol)
    {
      const CFreal etaSol = m_solPnts1D[iSol];
      m_etaFac[iSol] = 1.;
      for (CFuint iFac = 0; iFac < 1; ++iFac)
      {
        if (iFac != iSol)
        {
          const CFreal etaFac = m_solPnts1D[iFac];
          m_etaFac[iSol] *= (eta-etaFac)/(etaSol-etaFac);
        }
      }
    }

    // zta factors
    for (CFuint iSol = 0; iSol < 1; ++iSol)
    {
      const CFreal ztaSol = m_solPnts1D[iSol];
      m_ztaFac[iSol] = 1.;
      for (CFuint iFac = 0; iFac < 1; ++iFac)
      {
        if (iFac != iSol)
        {
          const CFreal ztaFac = m_solPnts1D[iFac];
          m_ztaFac[iSol] *= (zta-ztaFac)/(ztaSol-ztaFac);
        }
      }
    }

    // compute shapefunctions
    CFuint iFunc = 0;
    for (CFuint iKsi = 0; iKsi < 1; ++iKsi)
    {
      const CFreal ksiFac = m_ksiFac[iKsi];
      for (CFuint iEta = 0; iEta < 1; ++iEta)
      {
        const CFreal etaFac = m_etaFac[iEta];
        for (CFuint iZta = 0; iZta < 1; ++iZta, ++iFunc)
        {
          shapeFunc[iFunc] = ksiFac*etaFac*m_ztaFac[iZta];
        }
      }
    }
  }

   /**
   * Compute the Gradient of the base Function
   */
  static void computeGradientStates(
         const std::vector<RealMatrix>& jacob,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& grad)
  {
    throw Common::NotImplementedException (FromHere(),"The gradient of the FR base functions is not implemented (and should not be necessary...)");
  }

  /**
   * Computes the normal to a face at the given mapped coordinates,
   * scaled with the 'face Jacobian determinant'.
   * (Normal has the dimensionality of the Face + 1)
   */
  static void computeFaceJacobDetVectorAtMappedCoords(const std::vector<RealVector>& mappedCoord,
      const std::vector<Framework::Node*>& nodes,
      std::vector<RealVector>& normal)
  {
    throw Common::NotImplementedException (FromHere(),getName() + "::computeFaceJacobDetVectorAtMappedCoords()");
  }

  /**
   * Computes the normal to a given mapped coordinate plane, at the given mapped coordinates
   */
  static void computeMappedCoordPlaneNormal(const std::vector<CFuint>& planeIdx,
                                            const std::vector<RealVector>& mappedCoord,
                                            const std::vector<Framework::Node*>& nodes,
                                            std::vector<RealVector>& normal)
  {
    throw Common::NotImplementedException (FromHere(),getName() + "::computeMappedCoordPlaneNormal()");
  }

  static void computeJacobian(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& jacob)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

  static void computeJacobianPlus1D(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& jacob)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

  static void computeJacobianPlus2D(
         const std::vector<Framework::Node*>& nodes,
         const std::vector<RealVector>& mappedCoord,
               std::vector<RealMatrix>& jacob)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

  static void computeJacobianDeterminant(
         const std::vector<RealVector>& mappedCoord,
         const std::vector<Framework::Node*>& nodes,
               std::valarray<CFreal>& detJacobian)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

  static void computeJacobianDeterminantPlus1D(
         const std::vector<RealVector>& mappedCoord,
         const std::vector<Framework::Node*>& nodes,
               std::valarray<CFreal>& detJacobian)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

  static void computeJacobianDeterminantPlus2D(
         const std::vector<RealVector>& mappedCoord,
         const std::vector<Framework::Node*>& nodes,
               std::valarray<CFreal>& detJacobian)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

  static void computeFaceJacobianDeterminant(
          const std::vector<RealVector>& mappedCoord,
          const std::vector<Framework::Node*>& nodes,
          const Framework::IntegratorPattern& pattern,
                std::vector<RealVector>& faceJacobian);

  /**
   * Get the name of this base function
   */
  static const std::string getName()
  {
    return "FluxReconstructionPrismP0";
  }

  /**
   * Get the ID for the solution integration
   */
  static Framework::InterpolatorID getInterpolatorID()
  {
    return _interpolatorID;
  }

  /**
   * Get the ID for the solution integration
   */
  static void setInterpolatorID(const Framework::InterpolatorID& id)
  {
    _interpolatorID = id;
  }

  static CFreal computeVolume(const std::vector<Framework::Node*>& nodes)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

  static RealVector computeCentroid(const std::vector<Framework::Node*>& nodes)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

  /**
   * Compute base function in coordinate
   * @param coord contains the coordinates to be mapped
   * @param nodes contains the nodes
   * @return RealVector containing the Mapped Coordinates
   */
  static RealVector computeMappedCoordinates(const RealVector& coord, const std::vector<Framework::Node*>& nodes);

  static RealVector computeMappedCoordinatesPlus1D(const RealVector& coord, const std::vector<Framework::Node*>& nodes);

  static RealVector computeMappedCoordinatesPlus2D(const RealVector& coord, const std::vector<Framework::Node*>& nodes)
  {
    throw Common::NotImplementedException
      (FromHere(), getName()+"::computeMappedCoordinatesPlus2D()");
  }

  static std::vector<RealVector> computeAvgFaceNormals(const std::vector<Framework::Node*>& nodes)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

  static std::vector<RealVector> computeFaceNormals(const RealVector mappedCoord, const std::vector<Framework::Node*>& nodes)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

  static RealVector computeAvgCellNormal(const std::vector<Framework::Node*>& nodes)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

  static RealVector computeCellNormal(const RealVector& mappedCoord, const std::vector<Framework::Node*>& nodes)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

 /**
  * Check if a point (defined with mapped coordinates) is inside an element
  */
  static bool isInMappedElement(const RealVector& mappedCoord)
  {
    cf_assert(mappedCoord.size() == 3);
    if( (mappedCoord[0] >= -1.) &&
        (mappedCoord[0] <= +1.) &&
        (mappedCoord[1] >= -1.) &&
        (mappedCoord[1] <= +1.) &&
        (mappedCoord[2] >= -1.) &&
        (mappedCoord[2] <= +1.))
    {
      return true;
    }
    else
    {
      return false;
    }
  }

 /**
  * Check if a point is inside an element
  */
  static bool isInElement(const std::vector<Framework::Node*>& nodes, const RealVector& coord)
  {
    throw Common::NotImplementedException
      (FromHere(), getName()+"::isInElement()");
  }

private:

  /**
   * Default constructor without arguments
   */
  FluxReconstructionBaseFunctionPrismP0();

  /**
   * Default destructor
   */
  ~FluxReconstructionBaseFunctionPrismP0();

  static void computeFaceLineNormal(const CFuint& n,
                                    const std::vector<Framework::Node*>& nodes,
                                    const CFuint& i,
                                    const CFuint& j)
  {
    throw Common::ShouldNotBeHereException (FromHere(),"FR base functions should not be used as geometrical shape functions.");
  }

private: // data

  /// solution integrator ID
  static CFuint _interpolatorID;
  
  /// factors for computation of basis functions
  static RealVector m_ksiFac;

  /// factors for computation of basis functions
  static RealVector m_etaFac;

  /// factors for computation of basis functions
  static RealVector m_ztaFac;

  /// vector holding the 1D coordinates of solution points
  static RealVector m_solPnts1D;

}; // end of class FluxReconstructionBaseFunctionPrismP0

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_FluxReconstructionMethod_FluxReconstructionBaseFunctionPrismP0_hh
