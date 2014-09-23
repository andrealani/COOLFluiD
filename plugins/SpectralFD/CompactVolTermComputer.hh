#ifndef COOLFluiD_Numerics_SpectralFD_CompactVolTermComputer_hh
#define COOLFluiD_Numerics_SpectralFD_CompactVolTermComputer_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpectralFD/BaseVolTermComputer.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace SpectralFD {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a basic strategy that computes the volume terms in a cell with a fully compact scheme
 *
 * @author Kris Van den Abeele
 */
class CompactVolTermComputer : public BaseVolTermComputer {

public:  // types

  typedef Framework::BaseMethodStrategyProvider<
      SpectralFDMethodData,CompactVolTermComputer > PROVIDER;

public:  // methods

  /// Constructor
  CompactVolTermComputer(const std::string& name);

  /// Destructor
  ~CompactVolTermComputer();

  /// Gets the Class name
  static std::string getClassName()
  {
    return "CompactVolTermComputer";
  }

  /**
   * Set up private data and data
   */
  virtual void setup();

  /**
   * backup and reconstruct physical variable gradient in the required points
   */
  void backupAndReconstructGradients(const CFuint iVar, const std::vector< std::vector< RealVector >* >& cellGradients);

  /**
   * backup physical variable gradient in one cell in the required points
   */
  void backupPhysVarGrad(const CFuint iVar);

  /**
   * restore physical variable gradient in one cell in the required points
   */
  void restorePhysVarGrad(const CFuint iVar);

  /**
   * volume term contribution to a physical variable gradient
   */
  void computePhysVarGradientVolumeTerm(const CFuint iVar, std::vector< RealVector >& physVarGradUpdates);

protected: // functions

  /**
   * compute volume term contribution to the gradient of the given variable from the solution in the flux points
   */
  void computePhysVarGradVolTermFromFlxPntSol(const CFuint iVar, std::vector< RealVector >& physVarGradUpdates);

protected: // data

  /// backup for state gradients
  std::vector< RealVector > m_backupPhysVarGrad;

}; // class CompactVolTermComputer

//////////////////////////////////////////////////////////////////////////////

  }  // namespace SpectralFD

}  // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif  // COOLFluiD_Numerics_SpectralFD_CompactVolTermComputer_hh

