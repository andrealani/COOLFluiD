#ifndef COOLFluiD_Numerics_FluctSplit_STM_SplitterSys_hh
#define COOLFluiD_Numerics_FluctSplit_STM_SplitterSys_hh

//////////////////////////////////////////////////////////////////////////////

#include "SpaceTime_Splitter.hh"
#include "MathTools/RealMatrix.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MathTools { class MatrixInverter; }



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a generic interface for RDS splitters for
 * system schemes
 *
 * @author Thomas Wuilbaut
 *
 */
class STM_SplitterSys : public SpaceTime_Splitter {

public:

  /**
   * Constructor
   * @see Splitter()
   */
  STM_SplitterSys(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~STM_SplitterSys();

  /**
   * Set up
   */
  virtual void setup();

  /**
   * Compute the inflow parameters
   * @post K+ and K- will be computed
   */
  virtual void computeK(const std::vector<Framework::State*>& states,
			const InwardNormalsData* const normalsData);

private:
  
  /**
   * Sets the correct block limits for a System Splitter
   * Called by the STM_Splitter constructor.
   *
   * @see _nbEquations
   * @see _firstVarID
   * @see _lastVarID
   */
  void setBlockData();

protected: //data

  /// temporary data for holding identity matrix
  RealVector                       _identity;

  /// temporary data for holding eignvalues
  RealVector                       _eValues;

  /// temporary data for holding eignvaluesMinus
  RealVector                       _eValuesP;

  /// temporary data for holding eignvaluesPlus
  RealVector                       _eValuesM;

  /// temporary data for holding positive upwind parameter
  RealMatrix                       _tempKp;

  /// temporary data for holding negative upwind parameter
  RealMatrix                       _tempKm;

  /// temporary data for holding the matrix inverter
  MathTools::MatrixInverter*       _inverter;

  /// positive upwind parameters
  std::vector<RealMatrix*>         _kPlus;

  /// negative upwind parameters
  std::vector<RealMatrix*>         _kMin;

}; // end of class STM_SplitterSys

//////////////////////////////////////////////////////////////////////////////

    } // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_STM_SplitterSys_hh
