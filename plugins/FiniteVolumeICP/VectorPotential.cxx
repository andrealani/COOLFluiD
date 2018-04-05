#include "Common/NonCopyable.hh"
#include "MathTools/RealVector.hh"
#include "Framework/Framework.hh"
#include "FiniteVolumeICP/VectorPotential.hh"
#include "MathTools/MathConsts.hh"

/////////////////////////////////////////////////////////////////////////////

/*
using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::NavierStokes;
using namespace COOLFluiD::Physics::ICP;
using namespace COOLFluiD::Numerics::FiniteVolume;
using namespace COOLFluiD::MathTools;  
*/

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolumeICP {

//////////////////////////////////////////////////////////////////////////////

CFreal VectorPotential::getVectorPotentialRe()
{
    // don't need both... but what else can I do?
    CFreal vectorPotential_Re;
    CFreal vectorPotential_Im;

    // call function
    int vectorPotentialError = VectorPotentialFullParameters(vectorPotential_Re,vectorPotential_Im,
          _XtoCompute, _YtoCompute, _XelCurrent, _YelCurrent, 
          _elCurrent_Re, _elCurrent_Im, _permeability) ;
    if (vectorPotentialError==1) std::cout << "ERROR in Vector Potential computation\n";

    //return solution
    return vectorPotential_Re;
}

//////////////////////////////////////////////////////////////////////////////

CFreal VectorPotential::getVectorPotentialIm()
{
    // don't need both... but what else can I do?
    CFreal vectorPotential_Re;
    CFreal vectorPotential_Im;

    // call function
    int vectorPotentialError = VectorPotentialFullParameters(vectorPotential_Re,vectorPotential_Im,
          _XtoCompute, _YtoCompute, _XelCurrent, _YelCurrent, 
          _elCurrent_Re, _elCurrent_Im, _permeability) ;
    if (vectorPotentialError==1) std::cout << "ERROR in Vector Potential computation\n";

    //return solution
    return vectorPotential_Im;
}

//////////////////////////////////////////////////////////////////////////////

void VectorPotential::getVectorPotential(CFreal& vectorPotential_Re, CFreal& vectorPotential_Im) 
{
    // call function
    int vectorPotentialError = VectorPotentialFullParameters(vectorPotential_Re,vectorPotential_Im,
          _XtoCompute, _YtoCompute, _XelCurrent, _YelCurrent, 
          _elCurrent_Re, _elCurrent_Im, _permeability) ;
    if (vectorPotentialError==1) std::cout << "ERROR in Vector Potential computation\n";
}

//////////////////////////////////////////////////////////////////////////////

void VectorPotential::getElectricField(CFreal& electricField_Re, CFreal& electricField_Im) 
{
    CFreal vectorPotential_Re;
    CFreal vectorPotential_Im;

    // call function
    int vectorPotentialError = VectorPotentialFullParameters(vectorPotential_Re,vectorPotential_Im,
          _XtoCompute, _YtoCompute, _XelCurrent, _YelCurrent, 
          _elCurrent_Re, _elCurrent_Im, _permeability) ;
    if (vectorPotentialError==1) std::cout << "ERROR in Vector Potential computation\n";

    // real component of electric field intensity
    electricField_Re = 2.*MathTools::MathConsts::CFrealPi()*_frequency*vectorPotential_Im ;
    // real component of electric field intensity
    electricField_Im = -2.*MathTools::MathConsts::CFrealPi()*_frequency*vectorPotential_Re ;
}

//////////////////////////////////////////////////////////////////////////////

void VectorPotential::getElectricFieldStandardData(CFreal& electricField_Re, CFreal& electricField_Im) 
{
    CFreal vectorPotential_Re;
    CFreal vectorPotential_Im;

    // call function
    int vectorPotentialError = VectorPotentialFullParameters(vectorPotential_Re,vectorPotential_Im,
          _XtoCompute, _YtoCompute, _stdXelCurrent, _stdYelCurrent, 
          _stdElCurrent_Re, _stdElCurrent_Im, _stdPermeability) ;
    if (vectorPotentialError==1) std::cout << "ERROR in Vector Potential computation\n";

    // real component of electric field intensity
    electricField_Re = 2.*MathTools::MathConsts::CFrealPi()*_stdFrequency*vectorPotential_Im ;
    // real component of electric field intensity
    electricField_Im = -2.*MathTools::MathConsts::CFrealPi()*_stdFrequency*vectorPotential_Re ;
}

//////////////////////////////////////////////////////////////////////////////

int VectorPotential::VectorPotentialFullParameters(CFreal& vectorPotentialRe, CFreal& vectorPotentialIm,
 const CFreal& z, const CFreal& r,
 const std::vector<CFreal>& coilsR, const std::vector<CFreal>& coilsZ,
 const CFreal& elCurrentRe, const CFreal& elCurrentIm, const CFreal permeability) 
{
    CFreal zCoil=0.;
    CFreal rCoil=0.;
    CFreal k=0.;
    CFreal vectorPotential = 0.;
    CFuint nbCoils = coilsR.size();

    // some ghost cells could be below r=0....
    // let's apply odd symmetry about the axis!
    // by the way, we never apply 2D condition
    // to symmetry axis, so it's needed only computing LF.
   // CFreal radius = std::fabs(r);
     const CFreal RMIN = 1e-12;
     const CFreal rad = std::fabs(r);
     const CFreal radius = std::max(rad,RMIN);
     CFint rSign = r / radius;
        
      /*  if ( rSign < 0 ) {
            radius = r;
            rSign = +1;
      } */
    // prepare to & run elliptic integral:
    for (CFuint iCoil = 0; iCoil < nbCoils; ++iCoil) {
        rCoil = coilsR[iCoil];
        zCoil = coilsZ[iCoil];
        k = sqrt(4.*rCoil*radius/((rCoil+radius)*(rCoil+radius)+(z-zCoil)*(z-zCoil)));

        vectorPotential += sqrt(rCoil/std::max(radius,RMIN))*ellipticIntegralCombined(k);
    }

    // we need the real & imaginary part:
    vectorPotentialRe = rSign * vectorPotential *   .5*permeability*elCurrentRe/MathTools::MathConsts::CFrealPi();
    vectorPotentialIm = rSign * vectorPotential *   .5*permeability*elCurrentIm/MathTools::MathConsts::CFrealPi();

    // we assume that the simulation will be always axis-symmetric (r==0 is the axis)
    if (r==0.) {
        vectorPotentialRe = 0.; vectorPotentialIm = 0.;
    }

    // return error.
    if (std::isnan(vectorPotentialRe) || std::isnan(vectorPotentialIm))
       return 1;
    else
       return 0;
}

//////////////////////////////////////////////////////////////////////////////

CFreal VectorPotential::ellipticIntegralFirstKind(CFreal const& k)
{
  CFreal elIntFirstKind = 0.;
  const CFreal P = 1.- k*k;

  if(P < 1.e-15) elIntFirstKind = exp(1000.);
  else {
    elIntFirstKind = 1.38629436 + P*(0.096663443 + P*(0.035900924
         + P*(0.037425637 + 0.014511962*P)))
         - log(P)*(0.5 + P*(0.12498594 + P*(0.068802486
         + P*(0.033283553 + 0.0044178701*P))));
  }

  return elIntFirstKind;
}

//////////////////////////////////////////////////////////////////////////////

CFreal VectorPotential::ellipticIntegralSecondKind(CFreal const& k)
{
  CFreal elIntSecondKind = 0.;
  const CFreal P = 1.- k*k;

  if(P < 1.e-15) elIntSecondKind = 1.;
  else {
    elIntSecondKind = 1. + P*(0.44325141 + P*(0.062606012 + P*(0.047573836
         + 0.017365065*P)))
         - log(P)*P*(0.24998368 + P*(0.0920018
         + P*(0.040696975 + 0.0052644964*P)));
  }

  return elIntSecondKind;
}

//////////////////////////////////////////////////////////////////////////////

CFreal VectorPotential::ellipticIntegralCombined(CFreal const& k)
{
  CFreal elIntComb = 0.;

  if(k < 1.e-15) elIntComb = 0.;
  else elIntComb = ((2./k) - k)*ellipticIntegralFirstKind(k)
                 - (2./k)*ellipticIntegralSecondKind(k);

  return elIntComb;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

