#ifndef COOLFluiD_Numerics_FluctSplit_STKT_SplitterSyshyb_hh
#define COOLFluiD_Numerics_FluctSplit_STKT_SplitterSyshyb_hh

//////////////////////////////////////////////////////////////////////////////

#include "MathTools/RealMatrix.hh"
#include "FluctSplit/SpaceTime_Splitter.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools { class MatrixInverter; }



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a fluctuation splitter for system of equations
 * @todo NV:there is no moving mesh implemented there
 * @author Nadege Villedieu
 */
class STKT_SplitterSyshyb : public SpaceTime_Splitter {
public:

  /**
   * Constructor
   */
  STKT_SplitterSyshyb(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~STKT_SplitterSyshyb();

  /**
   * Set up
   */
  virtual void setup();

   /**
   * Compute the inflow parameters of a space-time element : k_space*dt/DIM + area/(DIM+1)
   * @post K+ and K- will be computed
   */
  virtual void computeK(const std::vector<Framework::State*>& states,
			const InwardNormalsData* const normalsData);

 
private:
  
  /**
   * Sets the correct block limits for a System Splitter
   * Called by the RDS_Splitter constructor.
   *
   * @see _nbEquations
   * @see _firstVarID
   * @see _lastVarID
   */
  void setBlockData();



protected: // data

  /// temporary data for holding identity matrix
  RealVector                       _identity;

  /// temporary data for holding eignvalues
  RealVector                       _eValues;

  /// temporary data for holding eignvaluesMinus
  RealVector                       _eValuesP;

  /// temporary data for holding eignvaluesPlus
  RealVector                       _eValuesM;

  /// temporary data for holding positive upwind parameter
  RealMatrix                       _rightEv;

  /// temporary data for holding negative upwind parameter
  RealMatrix                       _leftEv;

  /// temporary data for holding the matrix inverter
  MathTools::MatrixInverter*       _inverter;

  /// positive upwind parameters
  std::vector<RealMatrix*>         _kPlus;

  /// negative upwind parameters
  std::vector<RealMatrix*>         _kMin;

   /// Dimension
  CFuint DIM;

}; // end of class SpaceTimeRDS_SplitterSyshyb

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_STKT_SplitterSys_hh
