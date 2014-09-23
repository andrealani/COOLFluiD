#ifndef COOLFluiD_Numerics_FiniteVolume_FVMCC_Data_hh
#define COOLFluiD_Numerics_FiniteVolume_FVMCC_Data_hh

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {
      
//////////////////////////////////////////////////////////////////////////////

/**
 * This class represents a tuple grouping useful data for the FVM algorithm
 *
 * @author Andrea Lani
 *
 */
template <DeviceType DT, int ND, int NE>
class FVMCC_Data {
public:
  enum {DIM=ND, NBEQS=NE};
  
  typename MathTypes<DT>::VEC* stateL;
  
  typename MathTypes<DT>::VEC* stateR;
  
  typename MathTypes<DT>::VEC* nodeL;
  
  typename MathTypes<DT>::VEC* nodeR;
  
  typename MathTypes<DT>::VEC* n;
  
  typename MathTypes<DT, ND>::VEC unitNormal;
  
  typename MathTypes<DT, NE>::VEC flux;
  
  CFreal faceArea;
  
  CFuint stateIDL;
  
  CFuint stateIDR;
  
  CFreal* rhs;
  
  CFreal* updateCoeff;
  
  CFint isPerturb;
  
  CFint isBFace;
  
}; // end of class FVMCC_Data

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_Numerics_FiniteVolume_FVMCC_Data_hh
