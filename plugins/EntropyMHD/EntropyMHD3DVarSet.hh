#ifndef COOLFluiD_Physics_EntropyMHD_EntropyMHD3DVarSet_hh
#define COOLFluiD_Physics_EntropyMHD_EntropyMHD3DVarSet_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/ConvectiveVarSet.hh"
#include "Framework/EquationSetData.hh"
#include "EntropyMHDTerm.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace EntropyMHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * Base variable set for the EntropyMHD 3D entropy-variable MHD model with GLM.
 *
 * Eqs 0-6,8: standard Dedner et al. (2002) MHD + GLM.
 * Eq 7: entropy advection (flux = sigma*Vn, source = 0).
 * sigma = p*rho^(1-gamma); pressure recovered as p = sigma*rho^(gamma-1).
 *
 * Primitive: (rho, u, v, w, Bx, By, Bz, sigma, phi)
 *
 * @author Rayan Dhib
 */
class EntropyMHD3DVarSet : public Framework::ConvectiveVarSet {
public:

  typedef EntropyMHDTerm PTERM;

  EntropyMHD3DVarSet(Common::SafePtr<Framework::BaseTerm> term);

  virtual ~EntropyMHD3DVarSet();

  virtual void setup();

  virtual CFuint getBlockSeparator() const = 0;

  virtual void computeJacobians() = 0;

  /**
   * Entropy-variable MHD+GLM convective flux projected onto a normal direction.
   * Eq 7 = sigma*Vn (entropy advection). Fills _fluxArray (9 entries).
   */
  virtual void computeFlux(const RealVector& pdata,
                           const RealVector& normals);

  /**
   * Full 3D convective flux tensor (nEqs x nDim). Fills _physFlux.
   */
  virtual void computeStateFlux(const RealVector& pdata);

  /**
   * Compute the 9 eigenvalues of the GLM-MHD system.
   */
  virtual void computeEigenValues(const RealVector& pdata,
                                  const RealVector& normal,
                                  RealVector& eValues);

  virtual CFreal getMaxEigenValue(const RealVector& pdata,
                                  const RealVector& normal);

  virtual CFreal getMaxAbsEigenValue(const RealVector& pdata,
                                     const RealVector& normal);

  Common::SafePtr<EntropyMHDTerm> getModel() const
  {
    return _model;
  }

  virtual void setEigenVect1(RealMatrix& r,
                             RealMatrix& l,
                             RealVector& eValues,
                             const RealVector& normal) {}

  CFuint getNbEqs() const { return 9; }

  static std::vector<Framework::EquationSetData>& getEqSetData()
  {
    static std::vector<Framework::EquationSetData> eqSetData(1);
    return eqSetData;
  }

protected:

  /// acquaintance of the model
  Common::SafePtr<EntropyMHDTerm> _model;

}; // end of class EntropyMHD3DVarSet

//////////////////////////////////////////////////////////////////////////////

    } // namespace EntropyMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_EntropyMHD_EntropyMHD3DVarSet_hh
