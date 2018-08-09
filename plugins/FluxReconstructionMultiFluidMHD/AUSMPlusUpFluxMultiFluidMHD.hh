#ifndef COOLFluiD_FluxReconstructionMethod_AUSMPlusUpFluxMultiFluidMHD_hh
#define COOLFluiD_FluxReconstructionMethod_AUSMPlusUpFluxMultiFluidMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMultiFluidMHD/AUSMFluxMultiFluidMHD.hh"
#include "FluxReconstructionMultiFluidMHD/FluxReconstructionMultiFluidMHD.hh"
#include "FluxReconstructionMethod/FluxReconstructionSolverData.hh"
#include "FluxReconstructionMethod/RiemannFlux.hh"
#include "Framework/VarSetTransformer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace FluxReconstructionMethod {

//////////////////////////////////////////////////////////////////////////////

// This class computes the AUSM plus up flux

template <class UPDATEVAR>
class AUSMPlusUpFluxMultiFluid : public AUSMFluxMultiFluid<UPDATEVAR> {
public:

  /**
   * Constructor
   */
  AUSMPlusUpFluxMultiFluid(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~AUSMPlusUpFluxMultiFluid() ;
  
  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);

  /**
   * Set up private data
   */
  virtual void setup();

protected:

  /**
   * Compute the interface mass flux
   */
  virtual void computeMassFlux();

  /**
   * Compute the interface pressure flux
   */
  virtual void computePressureFlux();

  /// Correct the given Mach number
  virtual CFreal correctMachInf(CFreal oldMach) const
  {
    return oldMach;
  }

  virtual CFreal getCoeffKu() {return m_coeffKu;}
  virtual CFreal getCoeffKp() {return m_coeffKp;}
  virtual CFreal getCoeffSigma() {return m_coeffSigma;}
  virtual CFreal* getMachInf() {return &m_machInf[0];}
  virtual CFreal getBeta() {return m_beta;}
  virtual CFreal getFa() {return m_fa;}

private:

  /// preconditioning coefficient
  CFreal m_fa;

  /// user defined coefficient for Ku
  CFreal m_coeffKu;

  /// user defined coefficient for Kp
  CFreal m_coeffKp;

  /// user defined coefficient for sigma
  CFreal m_coeffSigma;

  /// mach infinity
  std::vector<CFreal> m_machInf;   //std::vector<CFreal>

  /// beta  coefficient
  CFreal m_beta;


}; // end of class AUSMPlusUpFlux

//////////////////////////////////////////////////////////////////////////////

  } // namespace FluxReconstructionMethod

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "FluxReconstructionMultiFluidMHD/AUSMPlusUpFluxMultiFluidMHD.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_AUSMPlusUpFluxMultiFluidMHD_hh
