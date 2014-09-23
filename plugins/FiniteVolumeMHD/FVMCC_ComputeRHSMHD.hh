#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRHSMHD_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRHSMHD_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/DataSocketSource.hh"
#include "Framework/DataSocketSink.hh"
#include "Common/CFMap.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {
    
    namespace FiniteVolume {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represent a command that computes the RHS and jacobian
 * using standard cell center FVM schemes for MHD simulations
 *
 * @author Andrea Lani
 * @author Mehmet Sarp Yalim
 *
 */
template <typename BASE>
class FVMCC_ComputeRHSMHD : public BASE {
public:

  /**
   * Defines the Config Option's of this class
   * @param options a OptionList where to add the Option's
   */
  static void defineConfigOptions(Config::OptionList& options);
    
  /**
   * Constructor.
   */
  explicit FVMCC_ComputeRHSMHD(const std::string& name);

  /**
   * Destructor.
   */
  virtual ~FVMCC_ComputeRHSMHD();

  /**
   * Configures this object by complementing the
   * implementation in ConfigObject
   */
  void configure ( Config::ConfigArgs& args );

  /**
   * Set up private data and data of the aggregated classes
   * in this command before processing phase
   */
  virtual void setup();

  /**
   * Unset up private data and data of the aggregated classes
   * in this command after processing phase
   */
  virtual void unsetup();

  /**
   * divBNodal values computed by weighted averaging
   */
  void computeDivBNodalValues(RealVector& divB);
  
  /**
   * Returns the DataSocket's that this command provides as sources
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSource> > providesSockets();

  /**
   * Returns the DataSocket's that this command needs as sinks
   * @return a vector of SafePtr with the DataSockets
   */
  std::vector<Common::SafePtr<Framework::BaseDataSocketSink> > needsSockets();
  
protected:
    
  /// Initialize the computation of RHS
  virtual void initializeComputationRHS();
  
  /// Compute between convective and diffusive term
  virtual void computeInterConvDiff();
  
  /// Finalize the computation of RHS
  virtual void finalizeComputationRHS();
  
private:

  /// socket for divB values at the nodes for outputting in Tecplot
  Framework::DataSocketSource<CFreal> socket_divBNodal;

  /// socket for divB values at the cell centers computed by Powell source term (if applicable)
  Framework::DataSocketSink<CFreal> socket_divBCellCenter;

  /// socket for average B values in x-direction at the cell faces to be used in the computation of divB in Powell source term (if applicable)
  Framework::DataSocketSource<CFreal> socket_avgBxFace;

  /// socket for average B values in y-direction at the cell faces to be used in the computation of divB in Powell source term (if applicable)
   Framework::DataSocketSource<CFreal> socket_avgByFace;
  
  /// socket for average B values in z-direction at the cell faces to be used in the computation of divB in Powell source term (if applicable)
  Framework::DataSocketSource<CFreal> socket_avgBzFace;
  
  /// flag checking of the first source term is "MHDConsACAST"
  bool _isFirstMHDConsACAST;
  
  /// name of the first source term computer
  std::string _firstSourceTermComputerName;
   
  /// Array storing the divB error values in the cell centers
  RealVector _divB;

  /// Output file save rate which should be EQUAL to the Tecplot file SaveRate
  CFuint _outputSaveRate;

  /// x coordinates of cell centers
  RealVector _currStateXCoord;
 
  /// y coordinates of cell centers
  RealVector _currStateYCoord;

  /// z coordinates of cell centers
  RealVector _currStateZCoord;

  /// corresponding to the vector of weights
  std::vector<CFreal>* _sumr;
                    
}; // class FVMCC_ComputeRHSMHD

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#include "FiniteVolumeMHD/FVMCC_ComputeRHSMHD.ci"

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_ComputeRhsJacobNumerics_hh
