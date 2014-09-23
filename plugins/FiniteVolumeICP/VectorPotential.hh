#ifndef COOLFluiD_Numerics_FiniteVolumeICP_VectorPotential_hh
#define COOLFluiD_Numerics_FiniteVolumeICP_VectorPotential_hh

//////////////////////////////////////////////////////////////////////////////

#include "Common/NonCopyable.hh"
#include "MathTools/RealVector.hh"
#include "Framework/Framework.hh"
#include <iostream>

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumeICP {

 /*   This class computes the Vector Potential
  *   usage:
  *    VectorPotential::getInstance().setPermeability(1000); 
  *
  *   @author Emanuele Sartori
  */
class VectorPotential : public Common::NonCopyable<VectorPotential> {

//////////////////////////////////////////////////////////////////////////////
//                      Class Vector Potential                              //
//                                                                          //
//  USAGE:                                                                  //
//  - set vacuum permeability and frequency of the current in coils         //
//  - set coils coordinates (can be a point or a vector of points)          //
//  - set coils current (can be only one value, Re and Im, not a vector)    //
//  - set coordinates of the point where we need the vector potential A     //
//  (you can set current and coil coords as standard)                       //
//                                                                          //
//  At this point you can use different methods:                            //
//  .getVectorPotentialRe();                                                //
//  .getVectorPotentialIm();                                                //
//  .getVectorPotential(vectorPotential_Re,vectorPotential_Im);             //
//  .getElectricField(electricField_Re,electricField_Im);                   //
//  .getElectricFieldStandardData(electricField_Re,electricField_Im);       //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

public:
  /// Default destructor
  ~VectorPotential()
  {
  }

  /// Returns the instance of this VectorPotential
  /// This is the access point to the Singleton
  /// @return the instance of the singleton
  static VectorPotential& getInstance()
  {
   static VectorPotential aVectorPotential;
   return aVectorPotential;
  }

  /*
   *  Set permeability
   */
  void setPermeability(CFreal permeability)
  {
    _permeability = permeability;
  }

  /*
   *  Set frequency
   */
  void setFrequency(CFreal frequency)
  {
    _frequency = frequency;
  }

  /*
   *  Set elCurrent, Re and Im components
   */
  void setElCurrent(CFreal elCurrent_Re, CFreal elCurrent_Im)
  {
    _elCurrent_Re = elCurrent_Re;
    _elCurrent_Im = elCurrent_Im;
  }

  /*
   * Set elCurrent coords, X and Y (z and r)
   * two versions:
   *  - setElCurrentCoords               if only one point has to be considered;
   *  - setElCurrentCoordsVectorCoords   if more points with the same el current 
   *                                     have to be considered.
   */
  void setElCurrentCoords(CFreal XelCurrent, CFreal YelCurrent)
  {
    _XelCurrent.resize(1);
    _YelCurrent.resize(1);
    _XelCurrent[0] = XelCurrent;
    _YelCurrent[0] = YelCurrent;
  }
  void promptElCurrentCoords()
  {
    std::cout << _XelCurrent[0] << "  ";
    std::cout << _YelCurrent[0] << "\n";
  }
  //  use this when you already have a vector of coordinates!!
  void setElCurrentVectorCoords(std::vector<CFreal> XelCurrent, std::vector<CFreal> YelCurrent)
  {
    // first resize _XelCurrent to XelCurrent.size()...
    _XelCurrent.resize( XelCurrent.size() );
    _YelCurrent.resize( YelCurrent.size() );
    _XelCurrent = XelCurrent;
    _YelCurrent = YelCurrent;
  }

  /*
   * coords where we need E: X and Y (z and r)
   */
  void setCoords(CFreal XtoCompute, CFreal YtoCompute)
  {
    _XtoCompute = XtoCompute;
    _YtoCompute = YtoCompute;
  }


  /*
   * Set standard coil data
   */
  void setStandardData()
  {
    _stdPermeability = _permeability;
    _stdFrequency    = _frequency;
    _stdElCurrent_Re = _elCurrent_Re;
    _stdElCurrent_Im = _elCurrent_Im;
    _stdXelCurrent.resize( _XelCurrent.size() );
    _stdYelCurrent.resize( _YelCurrent.size() );
    _stdXelCurrent = _XelCurrent;
    _stdYelCurrent = _YelCurrent;
  }

  /*
   * use standard coil data
   */
  void useStandardData()
  {
    _permeability = _stdPermeability;
    _frequency    = _stdFrequency;
    _elCurrent_Re = _stdElCurrent_Re;
    _elCurrent_Im = _stdElCurrent_Im;
    _XelCurrent.resize( _stdXelCurrent.size() );
    _YelCurrent.resize( _stdYelCurrent.size() );
    _XelCurrent = _stdXelCurrent;
    _YelCurrent = _stdYelCurrent;
  }

  /*
   * three different way to get Vector Potential !
   */
  CFreal getVectorPotentialRe();
  CFreal getVectorPotentialIm();
  void getVectorPotential(CFreal& vectorPotential_Re, CFreal& vectorPotential_Im);

  /*
   * and one to get E field
   */
  void getElectricField(CFreal& electricField_Re, CFreal& electricField_Im);

  /*
   * and one to get E field with only the coils contribution (standard data)
   */
  void getElectricFieldStandardData(CFreal& electricField_Re, CFreal& electricField_Im);

private:
  /// Constructor
  VectorPotential() {
  }

  /// let's compute Vector Potential 
  /**
   * This function returns Vector Potential (real or imaginary component) values
   * in some certain coordinates z,r
   * @param 
   */
  int VectorPotentialFullParameters(CFreal& vectorPotentialRe, CFreal& vectorPotentialIm,
     const CFreal& z, const CFreal& r,
     const std::vector<CFreal>& coilsR, const std::vector<CFreal>& coilsZ,
     const CFreal& elCurrentRe, const CFreal& elCurrentIm, const CFreal permeability) ;

  /**
   * This function returns values of the complete elliptic integral of the first kind
   * to the modulus k
   */
  CFreal ellipticIntegralFirstKind(CFreal const& k);
  /**
   * This function returns values of the complete elliptic integral of the second kind
   * to the modulus k
   */
  CFreal ellipticIntegralSecondKind(CFreal const& k);
  /**
   * This function returns values of the combination of the complete elliptic
   * integrals of the first and second kind to the modulus k
   */
  CFreal ellipticIntegralCombined(CFreal const& k);

private:

  // permeability
  CFreal _permeability;

  // frequency, in Hz
  CFreal _frequency;

  // coordinates where to compute E
  CFreal _XtoCompute;
  CFreal _YtoCompute;

  // coordinates where el current is running
  std::vector<CFreal> _XelCurrent;
  std::vector<CFreal> _YelCurrent;

  // Re and Im components of el current
  CFreal _elCurrent_Re;
  CFreal _elCurrent_Im;

  //// std data ////
  // permeability
  CFreal _stdPermeability;
  // frequency, in Hz
  CFreal _stdFrequency;
  // coordinates where el current is running
  std::vector<CFreal> _stdXelCurrent;
  std::vector<CFreal> _stdYelCurrent;
  // Re and Im components of el current
  CFreal _stdElCurrent_Re;
  CFreal _stdElCurrent_Im;
  //// ........ ////

}; // end class VectorPotential

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Framework_VectorPotential_hh
