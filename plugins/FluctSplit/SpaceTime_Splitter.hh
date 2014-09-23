#ifndef COOLFluiD_Numerics_FluctSplit_SpaceTime_Splitter_hh
#define COOLFluiD_Numerics_FluctSplit_SpaceTime_Splitter_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/SafePtr.hh"
#include "Common/NotImplementedException.hh"

#include "Framework/PhysicalModel.hh"
#include "Common/NullableObject.hh"
#include "Framework/ConvectiveVarSet.hh"

#include "FluctSplit/Splitter.hh"
#include "FluctSplit/InwardNormalsData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {



    namespace FluctSplit {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a generic RDS splitter
 *
 * @author Thomas Wuilbaut
 *
 */
class SpaceTime_Splitter : public Splitter {
public:

  /**
   * Default constructor without arguments
   */
  SpaceTime_Splitter(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~SpaceTime_Splitter();

  /**
   * Set up
   */
  virtual void setup();

  /**
   * Gets the Class name
   */
  static std::string getClassName()
  {
    return "SpaceTime_Splitter";
  }
  
  /**
   * Set the update coefficient
   */
  void setUpdateCoeff(Framework::DataHandle< CFreal> updateCoeff)
  {
    _updateCoeff = updateCoeff;
  }
  
  /**
   * Set the Volume of the cell at past and current states
   */
  void setCellVolume(const CFreal& p1)
  {
    _cellVolume = p1;
  }
  
  /**
   * Set the Volume of the cell at past and current states
   */
  void setPastCellVolume(const CFreal& p1)
  {
    _pastCellVolume = p1;
  }
  
  /**
   * Set the conservative States
   */
  void setConsStates(const std::vector<Framework::State*>& p1);

  /**
   * Set the intermediate States
   */
  void setInterStates(const std::vector<Framework::State*>& p1);

  /**
   * Set the intermediate conservative States
   */
  void setInterConsStates(const std::vector<Framework::State*>& p1);

  /**
   * Set the speed of the cell
   */
  void setCellSpeed(RealVector& p)
  {
    cf_assert(_cellSpeed.size() == p.size());
    _cellSpeed = p;
  }

  /**
   * Set the timestep
   */
  void setDT(const CFreal p){
    _timeStep = p;
  }

  /**
   * Distribute the residual due to the past
   */
  virtual void distributePast(const std::vector<Framework::State*>& tStates )
  {
     throw Common::NotImplementedException (FromHere(),"SpaceTime_Splitter::distributePast()");
  } 

 /**
   * Distribute the residual due to the past
   */
  virtual void distributePast( std::vector<RealVector>& residual )
  {
     throw Common::NotImplementedException (FromHere(),"SpaceTime_Splitter::distributePast()");
  } 

  /**
   * Compute the past dissipation of N scheme of Caraieni
   */
  virtual void  ComputePastDissipationAndTimeComp( const std::vector<Framework::State*>& tStates )
  {
     throw Common::NotImplementedException (FromHere(),"SpaceTime_Splitter::compute_past_dissipation()");
  }

  /**
   * Distribute the residual of intermediate states in K1
   */
  virtual void distributeInterK1(const std::vector<Framework::State*>& tStates,
                  std::vector<RealVector>& residual)
  {
     throw Common::NotImplementedException (FromHere(),"SpaceTime_Splitter::distributeInterK1()");
  }

  /**
   * Distribute the residual of intermediate states in K2
   */
  virtual void distributeInterK2(const std::vector<Framework::State*>& tStates,
                  std::vector<RealVector>& residual)
  {
     throw Common::NotImplementedException (FromHere(),"SpaceTime_Splitter::distributeInterK2()");
  }

protected: // data

  /// Store the Volume of the cell
  CFreal _cellVolume;

  /// Store the Past Volume of the cell
  CFreal _pastCellVolume;

  ///
  RealVector _nodeArea;

  /// Store the Conservative States
  std::vector<Framework::State*> _consStates;

  /// Store the Intermediate States
  std::vector<Framework::State*> _interStates;

  /// Store the Intermediate States
  std::vector<Framework::State*> _interConsStates;

  /// RealVector to store the cell speed
  RealVector _cellSpeed;

  /// Real to store the timeStep
  CFreal _timeStep;
  
  /// handle to current update coefficient storage
  Framework::DataHandle< CFreal> _updateCoeff;
  
   /// Pointer to the past residuals, contributed by the past states, computed with 1st order scheme
  /// RealVector is sized number of equations times number of states.
 std::vector<RealVector> past_residuals_order1;

}; // end of class SpaceTime_Splitter

//////////////////////////////////////////////////////////////////////////////

    }  // namespace FluctSplit



} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FluctSplit_SpaceTime_Splitter_hh
