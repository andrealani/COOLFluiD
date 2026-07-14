#ifndef COOLFluiD_Physics_EntropyMHD_EntropyMHD3DPrim_hh
#define COOLFluiD_Physics_EntropyMHD_EntropyMHD3DPrim_hh

//////////////////////////////////////////////////////////////////////////////

#include "EntropyMHD3DVarSet.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace EntropyMHD {

//////////////////////////////////////////////////////////////////////////////

/**
 * Primitive variable set for EntropyMHD 3D MHD with GLM projection.
 *
 * State: (rho, u, v, w, Bx, By, Bz, sigma, phi) where sigma = p*rho^(1-gamma).
 *
 * Eigensystem (computeEigenValuesVectors): the 9-wave Powell-Dedner GLM
 * eigenstructure is shared with standard MHD primitive (Roe 1981, Powell
 * 1994, Dedner 2002), but slot 7 holds the entropy variable sigma instead
 * of pressure p. The eigenvectors are obtained from the MHD primitive
 * eigenvectors by the closed-form basis transformation
 *
 *   T : (drho, dV, dB, dsigma, dphi) -> (drho, dV, dB, dp, dphi)
 *       T[p,rho]   = (gamma-1)*p/rho      (since p = sigma*rho^(gamma-1))
 *       T[p,sigma] = rho^(gamma-1)
 *       (identity elsewhere)
 *
 *   R_EntropyMHD = T^{-1} * R_MHD
 *   L_EntropyMHD = L_MHD * T
 *
 * which preserves the spectrum and only modifies row 7 of R and columns
 * 0,7 of L.
 *
 * @author Rayan Dhib
 */
class EntropyMHD3DPrim : public EntropyMHD3DVarSet {
public:

  EntropyMHD3DPrim(Common::SafePtr<Framework::BaseTerm> term);
  virtual ~EntropyMHD3DPrim();

  virtual void setup();

  CFuint getBlockSeparator() const { return 9; }
  void computeJacobians() {}

  void computePhysicalData(const Framework::State& state,
                           RealVector& data);

  void computeStateFromPhysicalData(const RealVector& data,
                                    Framework::State& state);

  /**
   * Right (R), left (L) eigenvector matrices and eigenvalue vector of the
   * EntropyMHD3D primitive Jacobian projected on the given normal.
   *
   * Wave ordering (EntropyMHD convention, matches EntropyMHD3DVarSet::computeEigenValues):
   *   k = 0 : Vn - cf   (fast magnetoacoustic, incoming if Vn < cf)
   *   k = 1 : Vn - ca   (Alfven, incoming if Vn < ca)
   *   k = 2 : Vn - cs   (slow magnetoacoustic, incoming if Vn < cs)
   *   k = 3 : Vn        (entropy / contact, sign of Vn)
   *   k = 4 : Vn + cs   (slow magnetoacoustic, outgoing if Vn > -cs)
   *   k = 5 : Vn + ca   (Alfven, outgoing if Vn > -ca)
   *   k = 6 : Vn + cf   (fast magnetoacoustic, outgoing if Vn > -cf)
   *   k = 7 : +c_h      (GLM cleaning, outgoing)
   *   k = 8 : -c_h      (GLM cleaning, incoming)
   *
   * Linearization is taken from the model's stored physical data
   * (getModel()->getPhysicalData()), as in MHD3DProjectionPrim.
   *
   * Implementation strategy: port MHD3DProjectionPrim::computeEigenValuesVectors
   * (Roe-Powell-Dedner GLM eigenvectors in primitive variables) and apply the
   * sigma-vs-p basis transformation (see class header). All entries except
   * row 7 of R and columns {0,7} of L are unchanged from the MHD form.
   */
  virtual void computeEigenValuesVectors(RealMatrix& rightEv,
                                         RealMatrix& leftEv,
                                         RealVector& eValues,
                                         const RealVector& normal);

  /**
   * Compute the eigensystem at an explicit physical-data linearization.
   * Same algorithm as computeEigenValuesVectors (which calls this with the
   * model's stored pdata), but accepts pdata as a parameter so callers can
   * compute eigenvectors at boundary states without polluting the model
   * linearization. This is the entry point for characteristic-based
   * boundary conditions that need the eigensystem at the boundary state.
   *
   * `linearData` must have layout matching EntropyMHDTerm enum (size getDataSize()).
   */
  void computeEigenValuesVectorsFromPData(const RealVector& linearData,
                                          RealMatrix& rightEv,
                                          RealMatrix& leftEv,
                                          RealVector& eValues,
                                          const RealVector& normal);

  void setDimensionalValues(const Framework::State& state,
                            RealVector& result)
  {
    result = static_cast<const RealVector&>(state);
  }

  void computePerturbedPhysicalData(const Framework::State& state,
                                    const RealVector& bData,
                                    RealVector& data,
                                    CFuint iVar)
  {
    computePhysicalData(state, data);
  }

  virtual void setStateVelocityIDs(std::vector<CFuint>& velIDs)
  {
    velIDs.resize(3);
    velIDs[XX] = 1;
    velIDs[YY] = 2;
    velIDs[ZZ] = 3;
  }

}; // end of class EntropyMHD3DPrim

//////////////////////////////////////////////////////////////////////////////

    } // namespace EntropyMHD

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Physics_EntropyMHD_EntropyMHD3DPrim_hh
