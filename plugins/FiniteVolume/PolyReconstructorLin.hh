#ifndef COOLFluiD_Numerics_FiniteVolume_PolyReconstructorLin_hh
#define COOLFluiD_Numerics_FiniteVolume_PolyReconstructorLin_hh

//////////////////////////////////////////////////////////////////////////////

#include "Framework/PolyReconstructor.hh"
#include "CellCenterFVMData.hh"

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

/**
 * This class offers a basic interface for polynomial reconstructors
 * when reconstructiong the states + the linearizedStates
 *
 * @author Thomas Wuilbaut
 *
 */
enum {LEFT=0, RIGHT=1};

class PolyReconstructorLin : public Framework::PolyReconstructor<CellCenterFVMData> {
public:

  /**
   * Constructor
   */
  PolyReconstructorLin(const std::string& name);

  /**
   * Default destructor
   */
  virtual ~PolyReconstructorLin();

  /**
   * Set private data that will be used during the computation
   */
  virtual void setup();

  /**
   * Configure the object
   */
  virtual void configure ( Config::ConfigArgs& args );

  /**
   * Set the pointers to extrapolated solution values
   */
  void setLeftValues2Ptr(std::vector<Framework::State*> *const values)
  {
    _extrapValues2[0] = values;
    _backupValues2[0].resize(values->size());
    for (CFuint i = 0; i < values->size(); ++i){
      _backupValues2[0][i].resize((*_extrapValues2[0])[i]->size());
    }
    cf_assert(_extrapValues2[0]->size() == _backupValues2[0].size());
    cf_assert(_extrapValues2[0] != CFNULL);
  }

  /**
   * Set the pointers to extrapolated solution values
   */
  void setRightValues2Ptr(std::vector<Framework::State*> *const values)
  {
    _extrapValues2[1] = values;
    _backupValues2[1].resize(values->size());
    for (CFuint i = 0; i < values->size(); ++i){
      _backupValues2[1][i].resize((*_extrapValues2[1])[i]->size());
    }
    cf_assert(_extrapValues2[1]->size() == _backupValues2[1].size());
    cf_assert(_extrapValues2[1] != CFNULL);
  }

  /**
   * Restore the back up values corresponding to the given variable ID
   */
  void restoreValues2(CFuint iVar, CFuint leftOrRight)
  {
    std::vector<Framework::State*>& values = getValues(leftOrRight);
    const std::vector<RealVector>& bValues = getBackupValues(leftOrRight);
    for (CFuint i = 0; i < bValues.size(); ++i) {
      (*values[i])[iVar] = bValues[i][iVar];
    }
    
    std::vector<Framework::State*>& values2 = getValues2(leftOrRight);
    const std::vector<RealVector>& bValues2 = getBackupValues2(leftOrRight);
    for (CFuint i = 0; i < bValues2.size(); ++i) {
      (*values2[i])[iVar] = bValues2[i][iVar];
    }
  }

 protected:

  /**
   * Get the left extrapolated solution values
   */
  std::vector<Framework::State*>& getValues2(CFuint leftOrRight)
  {
    cf_assert(_extrapValues2[leftOrRight] != CFNULL);
    return *_extrapValues2[leftOrRight];
  }

  /**
   * Get the back up extrapolated solution values
   */
  std::vector<RealVector>& getBackupValues2(CFuint leftOrRight)
  {
    return _backupValues2[leftOrRight];
  }


  /**
   * Copy the backup values LEFT and RIGHT
   */
  void copyBackupValues2(CFuint ip)
  {
    *getValues2(LEFT)[ip] = getBackupValues2(LEFT)[ip];
    *getValues2(RIGHT)[ip] = getBackupValues2(RIGHT)[ip];
  }
  
  /**
   * Extrapolate the solution in the face quadrature points
   */
  virtual void baseExtrapolateImpl(Framework::GeometricEntity* const face);
  
  /**
   * Extrapolate the solution in the face quadrature points
   */
  virtual void baseExtrapolateImpl(Framework::GeometricEntity* const face,
				   CFuint iVar, CFuint leftOrRight);
  
private:

  /// pointer to the extrapolated solution values
  std::vector<std::vector<Framework::State*>*> _extrapValues2;

  /// back up left values
  std::vector<std::vector<RealVector> > _backupValues2;

}; // end of class PolyReconstructorLin

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_PolyReconstructorLin_hh
