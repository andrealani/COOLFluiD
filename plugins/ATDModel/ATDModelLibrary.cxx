#include "ATDModel/ATDModelLibrary.hh"
#include "ATDModel/ATDModel.hh"
#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/BadValueException.hh"
#include "Framework/PhysicalModel.hh"
#include "Framework/MeshData.hh"
#include "Common/Stopwatch.hh"
#include "Environment/DirPaths.hh"
#include "Common/OSystem.hh"
#include "Common/StringOps.hh"

#include <fstream>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace ATDModel {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ATDModelLibrary,
                            PhysicalPropertyLibrary,
                            ATDModelModule,
                            1>
atdModelLibraryProvider("ATDModel");

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::defineConfigOptions(Config::OptionList& options)
{
  // dynamic options only for double, CFint, CFuint (unsigned integer)
  options.addConfigOption< CFdouble,
    Config::DynamicOption<> >("FactorOmega","Factor to reduce stiffness of chemical sorce terms.");

  options.addConfigOption< std::string >("mixtureName","Name of the mixture.");
  options.addConfigOption< std::string >("thermoName","Name of the thermodynamics' method.");
  options.addConfigOption< std::string >("chemName","Name of model for chemical reactions.");

  // to be removed
  options.addConfigOption< vector<CFuint> > ("ExampleIDs","IDs of variables or so.");

  options.addConfigOption< CFuint >
    ("DynamicViscosityModel","ID indentifying the viscosity model (1-Candler, 2-Yos/Gnoffo).");

  options.addConfigOption< CFuint >
    ("TempID","ID for Tq model (0-Tq=sqrt(T*Tv) 1-Tq=T^0.7*Tv^0.3).");

  options.addConfigOption< CFuint >
    ("RelaxID","Id of model for Qt-v: 0-(Zhong's Report) 1-(Millican and White) 2-(Zhong's report without using Tshk) 3-(Adapted from Mutation)");

  options.addConfigOption< CFuint >
    ("TetsID","Id of model for Ss: 0-Species dependent 1-All 5000");

  options.addConfigOption< CFdouble,
    Config::DynamicOption<>  >("Tshk","Temperature at shock.");

  options.addConfigOption< CFdouble >
    ("Tvsshk","Vibrational temperature at shock.");
}

//////////////////////////////////////////////////////////////////////////////

ATDModelLibrary::ATDModelLibrary(const std::string& name)
  : PhysicalChemicalLibrary(name),
    _ys(), _mmasses(), _hform(), _tau_vs(), _Tets(), _Ss(), _As(), _Bs(), _Cs(), _thermoCoefs(), _massTtmS(), _StchVecTemp(), _chemCoefsTemp(), _evibTemp(), _evibAsterTemp(), _omegaTemp(), _NIUsTemp(), _NIUsr(), _Asr()
{
  addConfigOptionsTo(this);
  
  _mixtureName = "";
  setParameter("mixtureName",&_mixtureName);

  _thermoName = "Candler";
  setParameter("thermoName",&_thermoName);

  _chemName = "air5park1";
  setParameter("chemName",&_chemName);

  _factorOmega = 1.0;
  setParameter("FactorOmega",&_factorOmega);

  _exampleIDs = vector<CFuint>();
  setParameter("ExampleIDs",&_exampleIDs);

  _etaModel = 1;
  setParameter("DynamicViscosityModel",&_etaModel);

  _TempID = 0;
  setParameter("TempID",&_TempID);

  _RelaxID = 0;
  setParameter("RelaxID",&_RelaxID);

  _TetsID = 0;
  setParameter("TetsID",&_TetsID);

  _Tshk = 14417.7;
  setParameter("Tshk",&_Tshk);

  _Tvsshk = 1833.0;
  setParameter("Tvsshk",&_Tvsshk);
}

//////////////////////////////////////////////////////////////////////////////

ATDModelLibrary::~ATDModelLibrary()
{
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::configure ( Config::ConfigArgs& args )
{
  PhysicalChemicalLibrary::configure(args);

  // here you are dataing in the input file settings
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::setup()
{
  PhysicalChemicalLibrary::setup();

  // if this is a parallel simulation, only ONE process at a time
  // sets the library
  setLibrary();
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::unsetup()
{
  if(isSetup()) {

    // deallocate local working arrays
    _ys.resize(0);
    _mmasses.resize(0);
    _hform.resize(0);
    _tau_vs.resize(0);
    _Tets.resize(0);
    _Ss.resize(0);
    _As.resize(0);
    _Bs.resize(0);
    _Cs.resize(0);
    _thermoCoefs.resize(0);
    _massTtmS.resize(0);
    _StchVecTemp.resize(0);
    _chemCoefsTemp.resize(0);
    _evibTemp.resize(0);
    _evibAsterTemp.resize(0);
    _CvvsTemp.resize(0);
    _CvvsAsterTemp.resize(0);
    _omegaTemp.resize(0);
    _NIUsTemp.resize(0);
    _NIUsr.resize(0,0);
    _Asr.resize(0,0);
    _mixtureName.resize(0);
    _thermoName.resize(0);
    _chemName.resize(0);
    _speciesNames.resize(0);
    _ChemReactArray.resize(0);
    _exampleIDs.resize(0);
    _molecIDs.resize(0);
    _CA.resize(0);
    _dissPartnersTemp.resize(0);
    _flagMoleculesIDs.resize(0);

    PhysicalChemicalLibrary::unsetup();
  }
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::getMolarMasses(RealVector& mm)
{
  // mm has size=nbSpecies
  mm=_mmasses;
}

//////////////////////////////////////////////////////////////////////////////

CFdouble ATDModelLibrary::getCvTr() const
{
  return _Rgas*_CvTrOR;
}

//////////////////////////////////////////////////////////////////////////////

CFdouble ATDModelLibrary::getMMass() const
{
    CFdouble massTtmTt=0.0;

    for (CFint i = 0; i < _NS; ++i) {
        massTtmTt+=_ys[i]/_mmasses[i];
    }

    return 1./massTtmTt;
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::setRiGas(RealVector& Ri)
{
  assert(Ri.size() == static_cast<CFuint>(_NS));
  for(CFint is = 0; is < _NS; ++is) {
    Ri[is] = _Rgas/_mmasses[is];
  }
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::setMoleculesIDs(std::vector<CFuint>& v)
{
  v.reserve(_NS);
  for (CFint i = 0; i < _NS; ++i) {
    if (_CA[i] > 1) {
      v.push_back(i);
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

CFdouble ATDModelLibrary::eta(CFdouble& temp,
                              CFdouble& pressure,
                              CFreal* tVec)
{
    CFdouble eta=0.0,PHIs;
    CFint i,j;

    for(i=0;i< _NS;i++)
    {
        _massTtmS[i]=_ys[i]/_mmasses[i];
        _NIUsTemp[i]=0.1*exp((_As[i]*log(temp)+_Bs[i])*log(temp)+_Cs[i]);
    }

    for(i=0;i< _NS;i++)
    {
        PHIs=0.0;
        for(j=0;j<_NS;j++)
        {
            PHIs+=_massTtmS[j]*pow(1+sqrt(_NIUsTemp[i]/_NIUsTemp[j])*pow(_mmasses[j]/_mmasses[i],0.25),2)/sqrt(8*(1+_mmasses[i]/_mmasses[j]));
        }

        eta+=_massTtmS[i]*_NIUsTemp[i]/PHIs;
    }

    return eta;
}

//////////////////////////////////////////////////////////////////////////////

CFdouble ATDModelLibrary::lambdaNEQ(CFdouble& temperature,
                                    CFdouble& pressure)
{
    CFint i,j;
    CFdouble PHIs,lambda=0.0;

    for(i=0;i< _NS;i++)
    {
        _massTtmS[i]=_ys[i]/_mmasses[i];
        _NIUsTemp[i]=0.1*exp((_As[i]*log(temperature)+_Bs[i])*log(temperature)+_Cs[i]);
    }

    for(i=0;i< _NS;i++)
    {
        PHIs=0.0;
        for(j=0;j<_NS;j++)
        {
            PHIs+=_massTtmS[j]*pow(1+sqrt(_NIUsTemp[i]/_NIUsTemp[j])*pow(_mmasses[j]/_mmasses[i],0.25),2)/sqrt(8*(1+_mmasses[i]/_mmasses[j]));
        }

        lambda+=_massTtmS[i]*_NIUsTemp[i]*(15.0/4+(_flagMoleculesIDs[i]? 1:0))/(_mmasses[i]*PHIs);
    }
    lambda*=_Rgas;

    return lambda;
}

//////////////////////////////////////////////////////////////////////////////

CFdouble ATDModelLibrary::lambdaEQ(CFdouble& temp, CFdouble& pressure)
{
    CFint i,j;
    CFdouble PHIs,lambda=0.0;

    for(i=0;i< _NS;i++)
    {
        _NIUsTemp[i]=0.1*exp((_As[i]*log(temp)+_Bs[i])*log(temp)+_Cs[i]);
    }

    for(i=0;i< _NS;i++)
    {
        PHIs=0.0;
        for(j=0;j<_NS;j++)
        {
            PHIs+=_massTtmS[j]*pow(1+sqrt(_NIUsTemp[i]/_NIUsTemp[j])*pow(_mmasses[j]/_mmasses[i],0.25),2)/sqrt(8*(1+_mmasses[i]/_mmasses[j]));
        }

        switch(_thermoID)
        {
        case 0:
            lambda+=_massTtmS[i]*_NIUsTemp[i]*(CvVibSpeciesOverR(temp,i)+(15.0/4+(_flagMoleculesIDs[i]? 1:0))/_mmasses[i])/PHIs;
            break;
        case 1:
            lambda+=_massTtmS[i]*_NIUsTemp[i]*(CpTeqOverRG(temp,i)+(5.0/4)/_mmasses[i])/PHIs;
            break;
        case 2:
            lambda+=_massTtmS[i]*_NIUsTemp[i]*(CpTeqOverRMB(temp,i)+(5.0/4)/_mmasses[i])/PHIs;
            break;
        }
    }
    lambda*=_Rgas;

    return lambda;
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::lambdaVibNEQ(CFreal& temperature,
                                   RealVector& tVec,
                                   CFdouble& pressure,
                                   CFreal& lambdaTrRo,
                                   RealVector& lambdaInt)
 {
    CFdouble PHIs;
    CFint i,j;

    lambdaTrRo=0.0;
    lambdaInt[0]=0.0;

    for(i=0;i< _NS;i++)
    {
        _massTtmS[i]=_ys[i]/_mmasses[i];
        _NIUsTemp[i]=0.1*exp((_As[i]*log(temperature)+_Bs[i])*log(temperature)+_Cs[i]);
    }

    for(i=0;i< _NS;i++)
    {
        PHIs=0.0;
        for(j=0;j<_NS;j++)
        {
            PHIs+=_massTtmS[j]*pow(1+sqrt(_NIUsTemp[i]/_NIUsTemp[j])*pow(_mmasses[j]/_mmasses[i],0.25),2)/sqrt(8*(1+_mmasses[i]/_mmasses[j]));
        }

        lambdaTrRo+=_massTtmS[i]*_NIUsTemp[i]*((5.0/2)*_atomicityCoeff[i]-3.0/2)/PHIs;

        switch(_thermoID)
        {
        case 0:
            lambdaInt[0]+=_massTtmS[i]*_NIUsTemp[i]*(CvVibSpeciesOverR(tVec[0],i))/PHIs;
            break;
        case 1:
            lambdaInt[0]+=_massTtmS[i]*_NIUsTemp[i]*(CpTeqOverRG(tVec[0],i)-(_atomicityCoeff[i]+1)/_mmasses[i])/PHIs;
            break;
        case 2:
            lambdaInt[0]+=_massTtmS[i]*_NIUsTemp[i]*(CpTeqOverRMB(tVec[0],i)-(_atomicityCoeff[i]+1)/_mmasses[i])/PHIs;
            break;
        }
    }

    lambdaTrRo*=_Rgas;
    lambdaInt[0]*=_Rgas;
}

//////////////////////////////////////////////////////////////////////////////

CFdouble ATDModelLibrary::sigma(CFdouble& temp,
                                CFdouble& pressure,
                                CFreal* tVec)
{
  /// AL: NOT to be implemented
  CFout << "Function not implemented: ATDModelLibrary::sigma() \n";
  throw Common::NotImplementedException(FromHere(),"ATDModelLibrary::sigma()");
  return 0;
}
CFdouble ATDModelLibrary::sigma_debug(CFdouble& temp, CFdouble& pressure, CFreal* tVec, CFuint elem_no)
{
  // Vatsalya
  CFout << "Function not implemented: ATDModelLibrary::sigma_debug() \n";
  throw Common::NotImplementedException(FromHere(),"ATDModelLibrary::sigma()");
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::gammaAndSoundSpeed(CFdouble& temp,
                                         CFdouble& pressure,
                                         CFdouble& rho,
                                         CFdouble& gamma,
                                         CFdouble& soundSpeed)
{
   /// AL: NOT to be implemented
  CFout << "Function not implemented: ATDModelLibrary::gammaAndSoundSpeed() \n";
  throw Common::NotImplementedException(FromHere(),"ATDModelLibrary::gammaAndSoundSpeed()");
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::frozenGammaAndSoundSpeed(CFdouble& temp,
                                                CFdouble& pressure,
                                                CFdouble& rho,
                                                CFdouble& gamma,
                                                CFdouble& soundSpeed,
                                                RealVector* tVec)
{
    _massTtmTt=0.0;
    _CvTrOR=0.0;

    for (CFint i = 0; i < _NS; ++i) {
        _massTtmTt+=_ys[i]/_mmasses[i];
        _CvTrOR+=_atomicityCoeff[i]*_massTtmS[i];
    }
    gamma = 1 + _massTtmTt/_CvTrOR;
    soundSpeed = sqrt(gamma*pressure/rho);
}

//////////////////////////////////////////////////////////////////////////////

CFdouble ATDModelLibrary::soundSpeed(CFdouble& temp, CFdouble& pressure)
{
  /// AL: NOT to be implemented
  CFout << "Function not implemented: ATDModelLibrary::soundSpeed() \n";
  throw Common::NotImplementedException(FromHere(),"ATDModelLibrary::soundSpeed()");
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::setComposition(CFdouble& temp,
                                     CFdouble& pressure,
                                     RealVector* x)
{
  /// AL: NOT to be implemented
  CFout << "Function not implemented: ATDModelLibrary::setComposition() \n";
  throw Common::NotImplementedException(FromHere(),"ATDModelLibrary::setComposition()");
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::resetComposition(const RealVector& x)
{ 
  /// AL: NOT to be implemented
  CFout << "Function not implemented: ATDModelLibrary::resetComposition() \n";
  throw Common::NotImplementedException(FromHere(),"ATDModelLibrary::resetComposition()");
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::setDensityEnthalpyEnergy(CFdouble& temp,
                                               CFdouble& pressure,
                                               RealVector& dhe)
{
  /// AL: NOT to be implemented
  CFout << "Function not implemented: ATDModelLibrary::setDensityEnthalpyEnergy() \n";
  throw Common::NotImplementedException(FromHere(),"ATDModelLibrary::setDensityEnthalpyEnergy()");
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::setDensityEnthalpyEnergy(CFdouble& temp,
                                                RealVector& tVec,
                                                CFdouble& pressure,
                                                RealVector& dhe,
                                                bool storeExtraData)
{
  CFdouble hformTot=0.0;
  CFint i;
    
  _massTtmTt=0.0;
  _CvTrOR=0.0;

  for ( i = 0; i < _NS; ++i) {
    _massTtmS[i]=_ys[i]/_mmasses[i];
    _massTtmTt+=_massTtmS[i];
    _CvTrOR+=_atomicityCoeff[i]*_massTtmS[i];
  }

  if (!storeExtraData) {
    dhe[0] = pressure/(_Rgas*_massTtmTt*temp);
    dhe[1] = 0.0;
    dhe[2] = 0.0;
    dhe[3] = 0.0;
    
    switch (_thermoID){

      case 0:

        for( i = 0; i < _NS; ++i) {
          dhe[3] += _ys[i]*energyVibSpeciesOverR(tVec[0],i);
          hformTot += _ys[i]*_hform[i];
        }
        dhe[2] = dhe[3]+_CvTrOR*temp;
        dhe[2] *= _Rgas;
        dhe[2] += hformTot;
        dhe[3] *= _Rgas;
        dhe[1] = dhe[2] + pressure/dhe[0];
	
        break;

      case 1:
        for( i = 0; i < _NS; ++i) {
          dhe[1] += _ys[i]*enthalpyTeqOverRG(tVec[0],i);
          hformTot += _ys[i]*_hform[i];
        }
        dhe[3] += dhe[1]-(_CvTrOR+_massTtmTt)*(tVec[0]-_Tref);
        dhe[3] *= _Rgas;
        dhe[3] -= hformTot;
        dhe[1] += (_CvTrOR+_massTtmTt)*(temp+_Tref-tVec[0]);

        dhe[1] *= _Rgas;

        dhe[2] = dhe[1] - pressure/dhe[0];

        break;

      case 2:

        for( i = 0; i < _NS; ++i) {
          dhe[1] += _ys[i]*enthalpyTeqOverRMB(tVec[0],i);
          hformTot += _ys[i]*_hform[i];
        }
        dhe[3] += dhe[1]-(_CvTrOR+_massTtmTt)*(tVec[0]-_Tref);
        dhe[3] *= _Rgas;
        dhe[3] -= hformTot;
        dhe[1] += (_CvTrOR+_massTtmTt)*(temp+_Tref-tVec[0]);

        dhe[1] *= _Rgas;

        dhe[2] = dhe[1] - pressure/dhe[0];

        break;
    }
  }
  else {
    CFdouble aux1,aux2;

    dhe[0] = pressure/(_Rgas*_massTtmTt*temp);
    dhe[1] = 0.0;
    dhe[2] = 0.0;
    dhe[3] = 0.0;

    switch (_thermoID){

      case 0:

        _extraData.cpVib=0.0;

        for(i = 0; i < _NS; ++i) {
          energyVibSpeciesOverR(tVec[0],i,aux1,aux2);
          dhe[3] += _ys[i]*aux1;
          hformTot += _ys[i]*_hform[i];

          _extraData.enthalpyTt[i]=_Rgas*((_atomicityCoeff[i]+1)*temp/_mmasses[i]+aux1)+_hform[i];
          _extraData.energyTr[i]=_Rgas*_atomicityCoeff[i]*temp/_mmasses[i]+_hform[i];
          _extraData.energyVib[i]=_Rgas*aux1;
          _extraData.enthalpyForm[i]=_hform[i];
          _extraData.dRhoEdRhoi[i]=_extraData.energyTr[i];
          _extraData.dRhoEvdRhoi[i]=_extraData.energyVib[i];
          _extraData.cpVib+=_ys[i]*aux2;
        }
        dhe[2] = dhe[3]+_CvTrOR*temp;
        dhe[2] *= _Rgas;
        dhe[2] += hformTot;
        dhe[3] *= _Rgas;
        dhe[1] = dhe[2] + pressure/dhe[0];

        _extraData.cpVib*=_Rgas;
        _extraData.dHdT=_Rgas*(_CvTrOR+_massTtmTt);
        _extraData.dEdT=_Rgas*_CvTrOR;

        break;

      case 1:

        _extraData.cpVib=0.0;

        for(i = 0; i < _NS; ++i) {
          enthalpyTeqOverRG(tVec[0],i,aux1,aux2);
          dhe[1] += _ys[i]*aux1;
          hformTot += _ys[i]*_hform[i];

          _extraData.enthalpyTt[i]=_Rgas*((_atomicityCoeff[i]+1)*(temp+_Tref-tVec[0])/_mmasses[i]+aux1);
          _extraData.energyTr[i]=_Rgas*_atomicityCoeff[i]*temp/_mmasses[i]+_hform[i];
          _extraData.energyVib[i]=_Rgas*(aux1-(_atomicityCoeff[i]+1)*(tVec[0]-_Tref)/_mmasses[i])-_hform[i];
          _extraData.enthalpyForm[i]=_hform[i];
          _extraData.dRhoEdRhoi[i]=_extraData.energyTr[i];
          _extraData.dRhoEvdRhoi[i]=_extraData.energyVib[i];
          _extraData.cpVib+=_ys[i]*(aux2-(_atomicityCoeff[i]+1)/_mmasses[i]);
        }
        dhe[3] += dhe[1]-(_CvTrOR+_massTtmTt)*(tVec[0]-_Tref);
        dhe[3] *= _Rgas;
        dhe[3] -= hformTot;
        dhe[1] += (_CvTrOR+_massTtmTt)*(temp+_Tref-tVec[0]);

        dhe[1] *= _Rgas;

        dhe[2] = dhe[1] - pressure/dhe[0];

        _extraData.cpVib*=_Rgas;
        _extraData.dHdT=_Rgas*(_CvTrOR+_massTtmTt);
        _extraData.dEdT=_Rgas*_CvTrOR;

        break;

      case 2:

        _extraData.cpVib=0.0;

        for(i = 0; i < _NS; ++i) {
          enthalpyTeqOverRMB(tVec[0],i,aux1,aux2);
          dhe[1] += _ys[i]*aux1;
          hformTot += _ys[i]*_hform[i];

          _extraData.enthalpyTt[i]=_Rgas*((_atomicityCoeff[i]+1)*(temp+_Tref-tVec[0])/_mmasses[i]+aux1);
          _extraData.energyTr[i]=_Rgas*_atomicityCoeff[i]*temp/_mmasses[i]+_hform[i];
          _extraData.energyVib[i]=_Rgas*(aux1-(_atomicityCoeff[i]+1)*(tVec[0]-_Tref)/_mmasses[i])-_hform[i];
          _extraData.enthalpyForm[i]=_hform[i];
          _extraData.dRhoEdRhoi[i]=_extraData.energyTr[i];
          _extraData.dRhoEvdRhoi[i]=_extraData.energyVib[i];
          _extraData.cpVib+=_ys[i]*(aux2-(_atomicityCoeff[i]+1)/_mmasses[i]);
        }
        dhe[3] += dhe[1]-(_CvTrOR+_massTtmTt)*(tVec[0]+_Tref-_Tref);
        dhe[3] *= _Rgas;
        dhe[3] -= hformTot;
        dhe[1] += (_CvTrOR+_massTtmTt)*(temp-tVec[0]);

        dhe[1] *= _Rgas;

        dhe[2] = dhe[1] - pressure/dhe[0];

        _extraData.cpVib*=_Rgas;
        _extraData.dHdT=_Rgas*(_CvTrOR+_massTtmTt);
        _extraData.dEdT=_Rgas*_CvTrOR;

        break;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

CFdouble ATDModelLibrary::density(CFdouble& temp,
                                   CFdouble& pressure,
                                   CFreal* tVec)
{
    _massTtmTt=0.0;

    for (CFint i = 0; i < _NS; ++i) {
        _massTtmTt+=_ys[i]/_mmasses[i];
    }
    return pressure/(_Rgas*_massTtmTt*temp);
}

//////////////////////////////////////////////////////////////////////////////

CFdouble ATDModelLibrary::energy(CFdouble& temp,
                                 CFdouble& pressure)
{
  /// AL: NOT to be implemented
  CFout << "Function not implemented: ATDModelLibrary::energy() \n";
  throw Common::NotImplementedException(FromHere(),"ATDModelLibrary::energy()");
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

CFdouble ATDModelLibrary::enthalpy(CFdouble& temp,
                                   CFdouble& pressure)
{
  /// AL: NOT to be implemented
    CFout << "Function not implemented: ATDModelLibrary::enthalpy() \n";
    throw Common::NotImplementedException(FromHere(),"ATDModelLibrary::setElemFractions()");
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::setElemFractions(const RealVector& yn)
{ 
  /// AL: NOT to be implemented
  CFout << "Function not implemented: ATDModelLibrary::setElemFractions() \n";
  throw Common::NotImplementedException(FromHere(),"ATDModelLibrary::setElemFractions()");
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::setElementXFromSpeciesY(const RealVector& ys)
{ 
  /// AL: NOT to be implemented
  CFout << "Function not implemented: ATDModelLibrary::setElementXFromSpeciesY() \n";
  throw Common::NotImplementedException(FromHere(),"ATDModelLibrary::setElementXFromSpeciesY()");
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::setElectronFraction(RealVector& ys)
{
  /// AL: NOT to be implemented
  CFout << "Function not implemented: ATDModelLibrary::setElectronFraction() \n";
  throw Common::NotImplementedException(FromHere(),"ATDModelLibrary::setElectronFraction()");
  // this affects cases with ionization
  // from Mutatiion
  //  // charge neutrality: xEl = sum(xIon)
  //   CFdouble yEl = 0.0;
  //   for (CFint is = 0; is < _NS; ++is) {
  //     if (_CHARGE[is] > 0) {
  //       yEl += ys[is] / _MOLARMASSP[is];
  //     }
  //   }

  //   yEl *= _MOLARMASSP[0]; // 1st species: electron
  //   ys[0] = yEl; // overwrite electron mass fraction
}

//////////////////////////////////////////////////////////////////////

void ATDModelLibrary::setSpeciesFractions(const RealVector& ys)
{
  for (CFint is = 0; is < _NS; ++is) {
    _ys[is] = ys[is];

    if (_ys[is] < 0.0) _ys[is] = 0.0;
    cf_assert(_ys[is] < 1.1);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::getSpeciesMolarFractions
(const RealVector& ys, RealVector& xs)
{
  CFdouble massTot=0.0;

  for (CFint is = 0; is < _NS; ++is) {
    cf_always_assert (ys[is] > 1.1);

    const CFreal mm = ys[is]/_mmasses[is];
    massTot += mm;
    xs[is] = mm;
  }
  xs *= 1./massTot;
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::getSpeciesMassFractions
(const RealVector& xs, RealVector& ys)
{
  CFdouble massTot=0.0;

  for (CFint is = 0; is < _NS; ++is) {
    cf_assert (xs[is] < 1.1);
    const CFreal mm = xs[is]*_mmasses[is];
    massTot += mm;
    ys[is] = mm;
  }
  ys /= massTot;
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::getSpeciesMassFractions(RealVector& ys)
{
  ys = _ys;
}

//////////////////////////////////////////////////////////////////////////////

CFdouble ATDModelLibrary::pressure(CFdouble& rho,
                                   CFdouble& temp,
                                   CFreal* tVec)
{
    _massTtmTt=0.0;

    for (CFint i = 0; i < _NS; ++i) {
        _massTtmTt+=_ys[i]/_mmasses[i];
    }
    return rho*_Rgas*_massTtmTt*temp;
}

//////////////////////////////////////////////////////////////////////////////

CFdouble ATDModelLibrary::electronPressure(CFreal rhoE,
                                           CFreal tempE)
{
  CFout << "Function not implemented: ATDModelLibrary::electronPressure() \n";
  throw Common::NotImplementedException(FromHere(),"ATDModelLibrary::electronPressure()");
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::getMassProductionTerm(CFdouble& temperature,
                                            RealVector& tVec,
                                            CFdouble& pressure,
                                            CFdouble& rho,
                                            const RealVector& ys,
                                            bool flagJac,
                                            RealVector& omega,
                                            RealMatrix& jacobian)
  {
    CFint i;
    omega=0.0;
    
    for(i=0;i<(CFint) _ChemReactArray.size();i++)
    {
        _ChemReactArray[i].omegaContribution(temperature,tVec,pressure,rho,ys,_mmasses,omega);
    }
    omega*=_mmasses;

    if (flagJac)
    {
        CFout << "Function not implemented: ATDModelLibrary::getMassProductionTerm(...,flagJac=true ,...) \n";
        throw Common::NotImplementedException(FromHere(),"ATDModelLibrary::getMassProductionTerm()");
    }
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::getSource(CFdouble& temperature,
                                 RealVector& tVec,
                                 CFdouble& pressure,
                                 CFdouble& rho,
                                 const RealVector& ys,
                                 bool flagJac,
                                 RealVector& omega,
                                 RealVector& omegav,
                                 CFdouble& omegaRad,
                                 RealMatrix& jacobian)
{
    CFint i;
    CFdouble aux1,aux2;

    if (flagJac)
    {
        CFout << "Function not implemented: ATDModelLibrary::getSource(...,flagJac=true ,...) \n";
        throw Common::NotImplementedException(FromHere(),"ATDModelLibrary::getSource()");
    }

    _massTtmTt=0.0;
    _CvTrOR=0.0;

    for ( i = 0; i < _NS; ++i) {
      cf_assert(ys[i] < 1.1);

      _massTtmS[i]=ys[i]/_mmasses[i];
      _massTtmTt+=_massTtmS[i];
      _CvTrOR+=_atomicityCoeff[i]*_massTtmS[i];
    }

    switch (_thermoID){

      case 0:

        for(i = 0; i < _NS; ++i) {
            energyVibSpeciesOverR(tVec[0],i,aux1,aux2);
            _evibTemp[i]=_Rgas*aux1;
            _CvvsTemp[i]=_Rgas*aux2;
            energyVibSpeciesOverR(temperature,i,aux1,aux2);
            _evibAsterTemp[i]=_Rgas*aux1;
            _CvvsAsterTemp[i]=_Rgas*aux2;
        }

        break;

      case 1:

        for(i = 0; i < _NS; ++i) {
            enthalpyTeqOverRG(tVec[0],i,aux1,aux2);
            _evibTemp[i]=_Rgas*(aux1-(_atomicityCoeff[i]+1)*(tVec[0]-_Tref)/_mmasses[i])-_hform[i];
            _CvvsTemp[i]=_Rgas*(aux2-(_atomicityCoeff[i]+1)/_mmasses[i]);
            enthalpyTeqOverRG(temperature,i,aux1,aux2);
            _evibTemp[i]=_Rgas*(aux1-(_atomicityCoeff[i]+1)*(temperature-_Tref)/_mmasses[i])-_hform[i];
            _CvvsAsterTemp[i]=_Rgas*(aux2-(_atomicityCoeff[i]+1)/_mmasses[i]);
        }

        break;

      case 2:

        for(i = 0; i < _NS; ++i) {
            enthalpyTeqOverRMB(tVec[0],i,aux1,aux2);
            _evibTemp[i]=_Rgas*(aux1-(_atomicityCoeff[i]+1)*(tVec[0]-_Tref)/_mmasses[i])-_hform[i];
            _CvvsTemp[i]=_Rgas*(aux2-(_atomicityCoeff[i]+1)/_mmasses[i]);
            enthalpyTeqOverRMB(temperature,i,aux1,aux2);
            _evibAsterTemp[i]=_Rgas*(aux1-(_atomicityCoeff[i]+1)*(temperature-_Tref)/_mmasses[i])-_hform[i];
            _CvvsAsterTemp[i]=_Rgas*(aux2-(_atomicityCoeff[i]+1)/_mmasses[i]);
        }

        break;
    }

    omegaRad=0.0;

    omega=0.0;

    for(i=0;i<(CFint) _ChemReactArray.size();i++)
    {
        _ChemReactArray[i].omegaContribution(temperature,tVec,pressure,rho,ys,_mmasses,omega);
    }
    omega*=_mmasses;

    ComputeOmegav(temperature,tVec,pressure,rho,ys,omega,omegav);
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::getRhoUdiff(CFdouble& temperature,
				  CFdouble& pressure,
				  RealVector& normConcGradients,
				  RealVector& normTempGradients,
				  CFreal* tVec,
				  RealVector& rhoUdiff,
				  bool fast)
{
  CFdouble PHIs,ktr=0.0,Cp=0.0;
  CFint i,j;

  for(i=0;i< _NS;i++)
  {

      _massTtmS[i]=_ys[i]/_mmasses[i];
      _NIUsTemp[i]=0.1*exp((_As[i]*log(temperature)+_Bs[i])*log(temperature)+_Cs[i]);
  }

  for(i=0;i< _NS;i++)
  {
      PHIs=0.0;
      for(j=0;j<_NS;j++)
      {
          PHIs+=_massTtmS[j]*pow(1+sqrt(_NIUsTemp[i]/_NIUsTemp[j])*pow(_mmasses[j]/_mmasses[i],0.25),2)/sqrt(8*(1+_mmasses[i]/_mmasses[j]));
      }

      ktr+=_massTtmS[i]*_NIUsTemp[i]*((5.0/2)*_atomicityCoeff[i]-3.0/2)/PHIs;

      switch(_thermoID)
      {
      case 0:
          Cp+=_ys[i]*((_atomicityCoeff[i]+1)/_mmasses[i]+CvVibSpeciesOverR(tVec[0],i));
          break;
      case 1:
          Cp+=_ys[i]*CpTeqOverRG(tVec[0],i);
          break;
      case 2:
          Cp+=_ys[i]*CpTeqOverRMB(tVec[0],i);
          break;
      }
  }
  Cp*=_Rgas;

  ktr*=_Rgas;
  for(i=0;i< _NS;i++)
  {   //Previuosly the formula was:
      //rhoUdiff[i]=-normConcGradients[i]*_Le*ktr/Cp;
      //but _Le should be in the denominator of the fraction
      rhoUdiff[i]=-normConcGradients[i]*ktr/(Cp*_Le);
  }
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::getDij_fick(RealVector& dx,
                                  CFdouble& pressure,
                                  CFdouble& temperature,
                                  RealMatrix& Dij,
                                  RealVector& rhoUdiff)
{
  //not implemented in mutation 2.0
  // look at plugins/MutationI/MutationLibrary.cxx
  CFout << "Function not implemented: ATDModelLibrary::getDij_fick() \n";
  throw Common::NotImplementedException(FromHere(),"ATDModelLibrary::getDij_fick()");
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::getGammaN(CFreal& m_GN)
{
  /// AL: NOT to be implemented
  CFout << "Function not implemented: ATDModelLibrary::getGammaN() \n";
  throw Common::NotImplementedException(FromHere(),"ATDModelLibrary::getGammaN()");
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::getGammaO(CFreal& m_GO)
{
  /// AL: NOT to be implemented
  CFout << "Function not implemented: ATDModelLibrary::getGammaO() \n";
  throw Common::NotImplementedException(FromHere(),"ATDModelLibrary::getGammaO()");
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::setSpeciesMolarFractions(const RealVector& xs)
{
  //not implemented in mutation 2.0
  CFout << "Function not implemented: ATDModelLibrary::setSpeciesMolarFractions() \n";
  throw Common::NotImplementedException(FromHere(),"ATDModelLibrary::setSpeciesMolarFractions()");
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::getSpeciesTotEnthalpies(CFdouble& temp,
                                               RealVector& tVec,
                                               CFdouble& pressure,
                                               RealVector& hsTot,
                                               RealVector* hsVib,
                                               RealVector* hsEl)
{
  CFint i, imol=0;
  CFdouble aux1;

  switch (_thermoID){

    case 0:
        for(i = 0; i < _NS; ++i) {
          aux1=energyVibSpeciesOverR(tVec[0],i);

          hsTot[i]=_Rgas*((_atomicityCoeff[i]+1)*temp/_mmasses[i]+aux1)+_hform[i];
          if(_flagMoleculesIDs[i])
          {
            (*hsVib)[imol]=_Rgas*aux1;
            imol++;
          }
        }

        break;

    case 1:
        for(CFint i = 0; i < _NS; ++i) {
          aux1=enthalpyTeqOverRG(tVec[0],i);

          hsTot[i]=_Rgas*((_atomicityCoeff[i]+1)*(temp+_Tref-tVec[0])/_mmasses[i]+aux1);
          if(_flagMoleculesIDs[i])
          {
            (*hsVib)[imol]=_Rgas*(aux1-(_atomicityCoeff[i]+1)*(tVec[0]-_Tref)/_mmasses[i])-_hform[i];
            imol++;
          }
        }

        break;

    case 2:
        for(CFint i = 0; i < _NS; ++i) {
          aux1=enthalpyTeqOverRMB(tVec[0],i);

          hsTot[i]=_Rgas*((_atomicityCoeff[i]+1)*(temp+_Tref-tVec[0])/_mmasses[i]+aux1);
          if(_flagMoleculesIDs[i])
          {
            (*hsVib)[imol]=_Rgas*(aux1-(_atomicityCoeff[i]+1)*(tVec[0]-_Tref)/_mmasses[i])-_hform[i];
            imol++;
          }
        }

        break;
  }
  
  
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::getSourceTermVT(CFdouble& temperature,
                                      RealVector& tVec,
                                      CFdouble& pressure,
                                      CFdouble& rho,
                                      RealVector& omegav,
                                      CFdouble& omegaRad)
{
    CFint i;
    CFdouble aux1,aux2;

    switch (_thermoID){

      case 0:

        for(i = 0; i < _NS; ++i) {
            energyVibSpeciesOverR(tVec[0],i,aux1,aux2);
            _evibTemp[i]=_Rgas*aux1;
            _CvvsTemp[i]=_Rgas*aux2;
            energyVibSpeciesOverR(temperature,i,aux1,aux2);
            _evibAsterTemp[i]=_Rgas*aux1;
            _CvvsAsterTemp[i]=_Rgas*aux2;
        }

        break;

      case 1:

        for(i = 0; i < _NS; ++i) {
            enthalpyTeqOverRG(tVec[0],i,aux1,aux2);
            _evibTemp[i]=_Rgas*(aux1-(_atomicityCoeff[i]+1)*(tVec[0]-_Tref)/_mmasses[i])-_hform[i];
            _CvvsTemp[i]=_Rgas*(aux2-(_atomicityCoeff[i]+1)/_mmasses[i]);
            enthalpyTeqOverRG(temperature,i,aux1,aux2);
            _evibTemp[i]=_Rgas*(aux1-(_atomicityCoeff[i]+1)*(temperature-_Tref)/_mmasses[i])-_hform[i];
            _CvvsAsterTemp[i]=_Rgas*(aux2-(_atomicityCoeff[i]+1)/_mmasses[i]);
        }

        break;

      case 2:

        for(i = 0; i < _NS; ++i) {
            enthalpyTeqOverRMB(tVec[0],i,aux1,aux2);
            _evibTemp[i]=_Rgas*(aux1-(_atomicityCoeff[i]+1)*(tVec[0]-_Tref)/_mmasses[i])-_hform[i];
            _CvvsTemp[i]=_Rgas*(aux2-(_atomicityCoeff[i]+1)/_mmasses[i]);
            enthalpyTeqOverRMB(temperature,i,aux1,aux2);
            _evibAsterTemp[i]=_Rgas*(aux1-(_atomicityCoeff[i]+1)*(temperature-_Tref)/_mmasses[i])-_hform[i];
            _CvvsAsterTemp[i]=_Rgas*(aux2-(_atomicityCoeff[i]+1)/_mmasses[i]);
        }

        break;
    }

    omegaRad=0.0;

    _omegaTemp=0.0;

    for(i=0;i<(CFint) _ChemReactArray.size();i++)
    {
        _ChemReactArray[i].omegaContribution(temperature,tVec,pressure,rho,_ys,_mmasses,_omegaTemp);
    }
    _omegaTemp*=_mmasses;

    ComputeOmegav(temperature,tVec,pressure,rho,_ys,_omegaTemp,omegav);
}

//////////////////////////////////////////////////////////////////////////////


void ATDModelLibrary::transportCoeffNEQ(CFreal& temperature,
					CFdouble& pressure,
					CFreal* tVec,
					RealVector& normConcGradients,
					RealVector& normTempGradients,
					CFreal& eta,
					CFreal& lambdaTrRo,
					RealVector& lambdaInt,
					RealVector& rhoUdiff)
{
    CFint i,j;
    CFdouble PHIs,Cp=0.0;

    lambdaTrRo=0.0;
    lambdaInt[0]=0.0;
    eta=0.0;

    for(i=0;i< _NS;i++)
    {

        _massTtmS[i]=_ys[i]/_mmasses[i];
        _NIUsTemp[i]=0.1*exp((_As[i]*log(temperature)+_Bs[i])*log(temperature)+_Cs[i]);
    }

    for(i=0;i< _NS;i++)
    {
        PHIs=0.0;
        for(j=0;j<_NS;j++)
        {
            PHIs+=_massTtmS[j]*pow(1+sqrt(_NIUsTemp[i]/_NIUsTemp[j])*pow(_mmasses[j]/_mmasses[i],0.25),2)/sqrt(8*(1+_mmasses[i]/_mmasses[j]));
        }

        lambdaTrRo+=_massTtmS[i]*_NIUsTemp[i]*(15.0/4+(_flagMoleculesIDs[i]? 1:0))/(_mmasses[i]*PHIs);

        eta+=_massTtmS[i]*_NIUsTemp[i]/PHIs;

        switch(_thermoID)
        {
        case 0:
            lambdaInt[0]+=_massTtmS[i]*_NIUsTemp[i]*(CvVibSpeciesOverR(tVec[0],i))/PHIs;
            Cp+=_ys[i]*((_atomicityCoeff[i]+1)/_mmasses[i]+CvVibSpeciesOverR(tVec[0],i));
            break;
        case 1:
            lambdaInt[0]+=_massTtmS[i]*_NIUsTemp[i]*(CpTeqOverRG(tVec[0],i)-(_atomicityCoeff[i]+1)/_mmasses[i])/PHIs;
            Cp+=_ys[i]*CpTeqOverRG(tVec[0],i);
            break;
        case 2:
            lambdaInt[0]+=_massTtmS[i]*_NIUsTemp[i]*(CpTeqOverRMB(tVec[0],i)-(_atomicityCoeff[i]+1)/_mmasses[i])/PHIs;
            Cp+=_ys[i]*CpTeqOverRMB(tVec[0],i);
            break;
        }
    }
    Cp*=_Rgas;

    lambdaTrRo*=_Rgas;
    lambdaInt[0]*=_Rgas;

    for(i=0;i< _NS;i++)
    {	//Previuosly the formula was:
        //rhoUdiff[i]=-normConcGradients[i]*_Le*lambdaTrRo/Cp;
        //but _Le should be in the denominator of the fraction
        rhoUdiff[i]=-normConcGradients[i]*lambdaTrRo/(Cp*_Le);
    }
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::getSourceEE(CFdouble& temperature,
                                   RealVector& tVec,
                                   CFdouble& pressure,
                                   CFdouble& rho,
                                   const RealVector& ys,
                                   bool flagJac,
                                   CFdouble& omegaEE)
{
  /// AL: NOT to be implemented
  CFout << "Function not implemented: ATDModelLibrary::getSourceEE() \n";
  throw Common::NotImplementedException(FromHere(),"ATDModelLibrary::getSourceEE()");
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::ReadDataMixture()
{
  string line;
  CFint linecount=0;

  _NS=0;

  ifstream mixfile((m_libPath+"mixtures/"+_mixtureName+".mix").c_str());

  if (mixfile.is_open())
  {

    while ( mixfile.good() )
    {
      getline(mixfile,line);

      // code to interpret the lines
      if (line[0] != '/')
      {
        linecount+=1;
        if (linecount==1)
        {
          _NS= static_cast<CFuint> (atoi(line.c_str()));

          // Arrays that will be reused often inside one or more of the functions and depend on _NS
          _ys.resize(_NS);
          _mmasses.resize(_NS);
          _hform.resize(_NS);
          _tau_vs.resize(_NS);
          _Tets.resize(_NS);
          _Ss.resize(_NS);
          _As.resize(_NS);
          _Bs.resize(_NS);
          _Cs.resize(_NS);
          _speciesNames.resize(_NS);
          _massTtmS.resize(_NS);
          _StchVecTemp.resize(_NS);
          _evibTemp.resize(_NS);
          _evibAsterTemp.resize(_NS);
          _CvvsTemp.resize(_NS);
          _CvvsAsterTemp.resize(_NS);
          _omegaTemp.resize(_NS);
          _NIUsTemp.resize(_NS);
          _NIUsr.resize(_NS,_NS);
          _Asr.resize(_NS,_NS);
          _CA.resize(_NS);

          _extraData.enthalpyTt.resize(_NS);
          _extraData.energyTr.resize(_NS);
          _extraData.energyVib.resize(_NS);
          _extraData.cpVib.resize(_nbTvib);
          _extraData.cpElec.resize(_NS);
          _extraData.eElec.resize(_NS);
          _extraData.enthalpyForm.resize(_NS);
          _extraData.dRhoEdRhoi.resize(_NS);
          _extraData.dRhoEvdRhoi.resize(_NS);
          _atomicityCoeff.resize(_NS);

          switch (_thermoID){
            case 0:
              _numTCR=3;
              _numTC=3;
              _thermoCoefs.resize(_numTC*_NS);
              break;
            case 1:
              _numTCR=6;
              _numTC=30;
              _thermoCoefs.resize(_numTC*_NS);
              break;
            case 2:
              _numTCR=8;
              _numTC=24;
              _thermoCoefs.resize(_numTC*_NS);
              break;
          }

          switch (_chemID){
            case 0:
              _chemCoefsTemp.resize(8);
              break;
            case 1:
              _chemCoefsTemp.resize(8);
              break;
            case 2:
              _chemCoefsTemp.resize(8);
              break;
            case 3:
              _chemCoefsTemp.resize(6);
              break;
          }

        }
        if (linecount>1 && linecount<_NS+2)  _speciesNames[linecount-2]=line;
      }
    }
    mixfile.close();
  }

  else {
    CFout << "Unable to open mixture file: " << m_libPath+"mixtures/"+_mixtureName+".mix" << "\n";
    abort();
  }
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::ReadDataChem()
{
  string line;
  ifstream chemfile((m_libPath+"chemistry/"+_chemName).c_str());
  bool getPartners=false;
  CFint idx;

  if (chemfile.is_open())
  {
    getline (chemfile,line);
    while ( strcoll(line.c_str(),"STOP") != 0 )
    {
      // code to interpret the lines
      if (line[0] != '/')
      {
          if(!getPartners){
              readReactSpecies(idx,_StchVecTemp,getPartners,line);
              readChemCoefs(idx,_chemCoefsTemp,line);

              if (!getPartners)
              {
                  ChemReact reaction(_StchVecTemp, _chemCoefsTemp, _chemID, _TempID);
                  _ChemReactArray.push_back(reaction);
              }
          }
          else{
              readPartners(_dissPartnersTemp,line);

              ChemReact reaction(_StchVecTemp, _chemCoefsTemp, _dissPartnersTemp, _chemID, _TempID);

              _ChemReactArray.push_back(reaction);

              getPartners=false;
          }
      }
      getline (chemfile,line);
    }
    chemfile.close();
  }
  else {
    CFout << "Unable to open chemistry file: " << m_libPath+"chemistry/"+_chemName << "\n";
    abort();
  }
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::ReadDataSpecies(CFint i)
{
  string line;
  CFint linecount=0;
  CFuint MMaux;
  stringstream StrStream;

  ifstream spcfile((m_libPath+"species/"+_speciesNames[i]+".spc").c_str());

  if (spcfile.is_open())
  {
    while ( spcfile.good() ){
      getline (spcfile,line);

      // code to interpret the lines
      if (line[0] != '/')
      {
        StrStream.str(line);
        linecount++;
        if ((linecount==1 && (StrStream>>_CA[i]).fail()) || (linecount==2 && (StrStream>>MMaux).fail())
            || (linecount==3 && (StrStream>>_hform[i]).fail()) || (linecount==4 && (StrStream>>_tau_vs[i]).fail())
            || (linecount==5 && (StrStream>>_Tets[i]).fail()) || (linecount==6 && (StrStream>>_As[i]).fail())
            || (linecount==7 && (StrStream>>_Bs[i]).fail()) || (linecount==8 && (StrStream>>_Cs[i]).fail()))
        {
          CFout << "Wrong format of species file: "<< m_libPath+"species/"+_speciesNames[i]+".spc" << "\n";
          abort();
        }

        StrStream.clear();
      }
    }
    _mmasses[i]=MMaux*1e-3;  // conversion from g to kg
    spcfile.close();
  }
  else {
    CFout << "Unable to open mixture file: " << m_libPath+"species/"+_speciesNames[i]+".spc" << "\n";
    abort();
  }

  if (linecount!=8)
  {
      CFout << "Wrong format of species file: " << m_libPath+"species/"+_speciesNames[i]+".spc" << "\n";
      abort();
  }
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::ReadDataThermo(CFint i)
{
  string line;
  int linecount=0;
  stringstream StrStream;

  ifstream thmfile((m_libPath+"thermo/"+_speciesNames[i]+_thermoName+".thm").c_str());

  if (thmfile.is_open())
  {
    CFint iind=i*_numTC;

    while ( thmfile.good() )
    {
      getline (thmfile,line);

      // code to interpret the lines
      if (line[0] != '/')
      {
        StrStream.str(line);
        linecount++;
        if (linecount>_numTC)
        {
            CFout << "Too many coefficients on file:"<< m_libPath+"thermo/"+_speciesNames[i]+_thermoName+".thm" << "\n";
            abort();
        }
        if ((StrStream>>_thermoCoefs[iind+linecount-1]).fail())
        {
          CFout << "Wrong format of thermo file:"<< m_libPath+"thermo/"+_speciesNames[i]+_thermoName+".thm" << "\n";
          abort();
        }
        StrStream.clear();
      }
    }
    thmfile.close();
  }
  else {
    CFout << "Unable to open mixture file" << "thermo/"+_speciesNames[i]+_thermoName+".thm" << "\n";
    abort();
  }

  if (linecount<_numTC)
  {
      CFout << "Not enough coefficients on file:"<< m_libPath+"thermo/"+_speciesNames[i]+_thermoName+".thm" << "\n";
      abort();
  }
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary:: setLibrary()
{
  
  CFint i,j;

  //universal constants
  _NA=6.0221418e23; //wikipedia
  _KB=1.3806503e-23; //google
  _PI=3.14159265;
  _Rgas=_NA*_KB;
  _Tref=298.16;
  _Le=1.4;

  //At the present stage we are only considering the two-temperature model and neutral species
  _nbTvib=1;
  _nbTe=0;

  CFout << "ATDModelLibrary::libpath     = "  << m_libPath << "\n";
  CFout << "ATDModelLibrary::mixtureName = "  << _mixtureName << "\n";
  CFout << "ATDModelLibrary::chemName = " << _chemName << "\n";
  CFout << "ATDModelLibrary::thermoName = " << _thermoName << "\n";

  // if nobody has configured the library path, it is done here
  // with default = basedir + "plugins/ATDModel/data/"
  if (m_libPath == "") {
    CFLog(NOTICE, "ATDModelLibrary::libpath set to default" << "\n");
    std::string baseDir = Environment::DirPaths::getInstance().getBaseDir().string();
    m_libPath = baseDir + "/plugins/ATDModel/data/";
  }

  // Define _thermoID
  if (_thermoName=="Candler") _thermoID=0;
  else if (_thermoName=="Gnoffo") _thermoID=1;
  else if (_thermoName=="McBride") _thermoID=2;
  else {
      CFout << "ThermoName" + _thermoName + "is not valid." << "\n";
      abort();
  }

  // Define _chemID
  if (_chemName=="air5park1" || _chemName=="nitrogen2park1") _chemID=0;
  else if (_chemName=="air5park2" || _chemName=="nitrogen2park2") _chemID=1;
  else if (_chemName=="air5park3" || _chemName=="nitrogen2park3") _chemID=2;
  else if (_chemName=="air5DunnKang" || _chemName=="nitrogen2DunnKang") _chemID=3;
  else {
      CFout << "chemName " + _chemName + " is not valid." << "\n";
      abort();
  }

  // Test _TempID
  if (_TempID!=0 && _TempID!=1) {
      CFout << "TempID " << _TempID << " is not valid." << "\n";
      abort();
  }

  // Test _RelaxID
  if (_RelaxID!=0 && _RelaxID!=1 && _RelaxID!=2 && _RelaxID!=3) {
      CFout << "RelaxID " << _RelaxID << " is not valid." << "\n";
      abort();
  }

  // Test _TetsID
  if (_TetsID!=0 && _TetsID!=1) {
      CFout << "TetsID " << _TetsID << " is not valid." << "\n";
      abort();
  }

  //Read data from files initialize arrays
  ReadDataMixture();

  ReadDataChem();

  CFout << "Mixture and chemical reactions' data read. \n";

  //Read data from files
  for (i = 0; i < _NS; ++i) {
    ReadDataSpecies(i);
    ReadDataThermo(i);
  }

  setMoleculesIDs(_molecIDs);

  _flagMoleculesIDs.resize(_NS,false);
  for (CFuint i = 0; i < _molecIDs.size(); ++i) {
    _flagMoleculesIDs[_molecIDs[i]] = true;
  }
  for (CFint i = 0; i < _NS; ++i) {
    _atomicityCoeff[i] = (_flagMoleculesIDs[i]) ? 2.5 : 1.5;
  }

  for(i=0;i<_NS;i++)
  {
      for(j=0;j<_NS;j++)
      {
          _NIUsr(i,j)=1000*(_mmasses[i]*_mmasses[j])/(_mmasses[i]+_mmasses[j]);
          _Asr(i,j)=(_tau_vs[i]!=-1)? 1.16e-3*sqrt(_NIUsr(i,j))*pow(_tau_vs[i],4.0/3):0;
      }
  }

  CFout << "Library set." << "\n";
  
}

//////////////////////////////////////////////////////////////////////////////

CFdouble ATDModelLibrary::energyVibSpeciesOverR(CFdouble& tVib,CFint& i)
{
  CFdouble engyOR=0.0;

  CFint iind=i*_numTC;

  if (_tau_vs[i] != -1) engyOR += _tau_vs[i]/(exp(_tau_vs[i]/tVib)-1);
  if (_thermoCoefs[iind] != -1) engyOR += _thermoCoefs[iind+2]*_thermoCoefs[iind+1]*exp(-_thermoCoefs[iind]/tVib)/(_thermoCoefs[iind+1]+_thermoCoefs[iind+2]*exp(-_thermoCoefs[iind]/tVib));

  return engyOR/_mmasses[i];
}

//////////////////////////////////////////////////////////////////////////////

CFdouble ATDModelLibrary::enthalpyTeqOverRG(CFdouble& temp,CFint& i)
{
  CFdouble etpyOR=0.0;

  CFint iind;

  if (temp<=300.){
//      CFout << "Warning: Temperature outside allowable range, " << temp << "K is less than 300K ( at enthalpyTeqOverRG() ). Extrapolating. \n";
      iind=i*_numTC;
  }
  else if (temp<=1000.) iind=i*_numTC;
  else if (temp<=6000.) iind=i*_numTC+_numTCR;
  else if (temp<=15000.) iind=i*_numTC+2*_numTCR;
  else if (temp<=25000.) iind=i*_numTC+3*_numTCR;
  else if (temp<=35000.) iind=i*_numTC+4*_numTCR;
  else {
//      CFout << "Warning: Temperature outside allowable range, " << temp << "K is greater than 35000K ( at enthalpyTeqOverRG() ). Extrapolating. \n";
      iind=i*i*_numTC+4*_numTCR;;
  }

  for (CFint j = 5; j >= 1; --j) {
    etpyOR += _thermoCoefs[iind+j-1]/j;
    etpyOR *= temp;
  }

  etpyOR += _thermoCoefs[iind+5];

  return etpyOR/_mmasses[i];
}

//////////////////////////////////////////////////////////////////////////////

CFdouble ATDModelLibrary::enthalpyTeqOverRMB(CFdouble& temp,CFint& i)
{
    CFdouble etpyOR=0.0;

    CFint iind;

    if (temp<=200.){
//        CFout << "Warning: Temperature outside allowable range, " << temp << "K is less than 200K ( at enthalpyTeqOverRMB() ). Extrapolating. \n";
        iind=i*_numTC;
    }
    else if (temp<=1000.) iind=i*_numTC;
    else if (temp<=6000.) iind=i*_numTC+_numTCR;
    else if (temp<=20000.) iind=i*_numTC+2*_numTCR;
    else {
//        CFout << "Warning: Temperature outside allowable range, " << temp << "K is greater than 20000K ( at enthalpyTeqOverRMB() ). Extrapolating. \n";
        iind=i*_numTC+2*_numTCR;
    }


    for (CFint j = 5; j >= 1; --j) {
      etpyOR += _thermoCoefs[iind+j+1]/j;
      etpyOR *= temp;
    }

    etpyOR += _thermoCoefs[iind+7];
    etpyOR -= _thermoCoefs[iind]/temp;
    etpyOR += _thermoCoefs[iind+1]*log(temp);

    return etpyOR/_mmasses[i];
}

//////////////////////////////////////////////////////////////////////////////

CFdouble ATDModelLibrary::CvVibSpeciesOverR(CFdouble& tVib,CFint& i)
{
    CFint iind=i*_numTC;
    CFdouble CvvOR=0.0;
    if (_tau_vs[i] != -1){
        CvvOR += pow(_tau_vs[i]/tVib,2)*exp(_tau_vs[i]/tVib)/pow(exp(_tau_vs[i]/tVib)-1,2);
    }
    if (_thermoCoefs[iind] != -1){
        CvvOR += _thermoCoefs[iind+2]*_thermoCoefs[iind+1]*pow(_thermoCoefs[iind]/tVib,2)*exp(-_thermoCoefs[iind]/tVib)/pow(_thermoCoefs[iind+1]+_thermoCoefs[iind+2]*exp(-_thermoCoefs[iind]/tVib),2);
    }

    return CvvOR/_mmasses[i];
}

//////////////////////////////////////////////////////////////////////////////

CFdouble ATDModelLibrary::CpTeqOverRG(CFdouble& temp,CFint& i)
{
    CFint iind;
    CFdouble CpOR=0.0;

    if (temp<=300.){
//        CFout << "Warning: Temperature outside allowable range, " << temp << "K is less than 300K ( at enthalpyTeqOverRG() ). Extrapolating. \n";
        iind=i*_numTC;
    }
    else if (temp<=1000.) iind=i*_numTC;
    else if (temp<=6000.) iind=i*_numTC+_numTCR;
    else if (temp<=15000.) iind=i*_numTC+2*_numTCR;
    else if (temp<=25000.) iind=i*_numTC+3*_numTCR;
    else if (temp<=35000.) iind=i*_numTC+4*_numTCR;
    else {
//        CFout << "Warning: Temperature outside allowable range, " << temp << "K is greater than 35000K ( at enthalpyTeqOverRG() ). Extrapolating. \n";
        iind=i*_numTC+4*_numTCR;
    }

    for (CFint j = 5; j >= 1; --j) {
      CpOR *= temp;
      CpOR += _thermoCoefs[iind+j-1];
    }

    return CpOR/_mmasses[i];
}

//////////////////////////////////////////////////////////////////////////////

CFdouble ATDModelLibrary::CpTeqOverRMB(CFdouble& temp,CFint& i)
{
    CFint iind;
    CFdouble CpOR=0.0;

    if (temp<=200.){
//        CFout << "Warning: Temperature outside allowable range, " << temp << "K is less than 200K ( at enthalpyTeqOverRMB() ). Extrapolating. \n";
        iind=i*_numTC;
    }
    else if (temp<=1000.) iind=i*_numTC;
    else if (temp<=6000.) iind=i*_numTC+_numTCR;
    else if (temp<=20000.) iind=i*_numTC+2*_numTCR;
    else {
//        CFout << "Warning: Temperature outside allowable range, " << temp << "K is greater than 20000K ( at enthalpyTeqOverRMB() ). Extrapolating. \n";
        iind=i*_numTC+2*_numTCR;
    }

    for (CFint j = 5; j >= 1; --j) {
      CpOR *= temp;
      CpOR += _thermoCoefs[iind+j+1];
    }

    CpOR += (_thermoCoefs[iind]/temp+_thermoCoefs[iind+1])/temp;

    return CpOR/_mmasses[i];
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::energyVibSpeciesOverR(CFdouble& tVib,CFint& i,CFdouble& EvOR,CFdouble &CvvOR)
{
  CFint iind=i*_numTC;

  EvOR=0.0;
  CvvOR=0.0;
  if (_tau_vs[i] != -1){
      EvOR += _tau_vs[i]/(exp(_tau_vs[i]/tVib)-1);
      CvvOR += pow(_tau_vs[i]/tVib,2)*exp(_tau_vs[i]/tVib)/pow(exp(_tau_vs[i]/tVib)-1,2);
  }
  if (_thermoCoefs[iind] != -1){
      EvOR += _thermoCoefs[iind+2]*_thermoCoefs[iind+1]*exp(-_thermoCoefs[iind]/tVib)/(_thermoCoefs[iind+1]+_thermoCoefs[iind+2]*exp(-_thermoCoefs[iind]/tVib));
      CvvOR += _thermoCoefs[iind+2]*_thermoCoefs[iind+1]*pow(_thermoCoefs[iind]/tVib,2)*exp(-_thermoCoefs[iind]/tVib)/pow(_thermoCoefs[iind+1]+_thermoCoefs[iind+2]*exp(-_thermoCoefs[iind]/tVib),2);
  }

  EvOR/=_mmasses[i];
  CvvOR/=_mmasses[i];
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::enthalpyTeqOverRG(CFdouble& temp,CFint& i,CFdouble& htotOR,CFdouble &CpOR)
{
  CFint iind;

  htotOR=0.0;
  CpOR=0.0;

  if (temp<=300.){
//      CFout << "Warning: Temperature outside allowable range, " << temp << "K is less than 300K ( at enthalpyTeqOverRG() ). Extrapolating. \n";
      iind=i*_numTC;
  }
  else if (temp<=1000.) iind=i*_numTC;
  else if (temp<=6000.) iind=i*_numTC+_numTCR;
  else if (temp<=15000.) iind=i*_numTC+2*_numTCR;
  else if (temp<=25000.) iind=i*_numTC+3*_numTCR;
  else if (temp<=35000.) iind=i*_numTC+4*_numTCR;
  else {
//      CFout << "Warning: Temperature outside allowable range, " << temp << "K is greater than 35000K ( at enthalpyTeqOverRG() ). Extrapolating. \n";
      iind=i*_numTC+4*_numTCR;
  }
  for (CFint j = 5; j >= 1; --j) {
    htotOR += _thermoCoefs[iind+j-1]/j;
    htotOR *= temp;
    CpOR *= temp;
    CpOR += _thermoCoefs[iind+j-1];
  }

  htotOR += _thermoCoefs[iind+5];

  CpOR/=_mmasses[i];
  htotOR/=_mmasses[i];
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::enthalpyTeqOverRMB(CFdouble& temp,CFint& i,CFdouble& htotOR,CFdouble &CpOR)
{
    CFint iind;

    htotOR=0.0;
    CpOR=0.0;

    if (temp<=200.){
//        CFout << "Warning: Temperature outside allowable range, " << temp << "K is less than 200K ( at enthalpyTeqOverRMB() ). Extrapolating. \n";
        iind=i*_numTC;
    }
    else if (temp<=1000.) iind=i*_numTC;
    else if (temp<=6000.) iind=i*_numTC+_numTCR;
    else if (temp<=20000.) iind=i*_numTC+2*_numTCR;
    else {
//        CFout << "Warning: Temperature outside allowable range, " << temp << "K is greater than 20000K ( at enthalpyTeqOverRMB() ). Extrapolating. \n";
        iind=i*_numTC+2*_numTCR;
    }


    for (CFint j = 5; j >= 1; --j) {
      htotOR += _thermoCoefs[iind+j+1]/j;
      htotOR *= temp;
      CpOR *= temp;
      CpOR += _thermoCoefs[iind+j+1];
    }

    htotOR += _thermoCoefs[iind+7]-_thermoCoefs[iind]/temp+_thermoCoefs[iind+1]*log(temp);
    CpOR += (_thermoCoefs[iind]/temp+_thermoCoefs[iind+1])/temp;

    htotOR/=_mmasses[i];
    CpOR/=_mmasses[i];
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::readReactSpecies(CFint& idx,RealVector& StchVec,bool& getPartners,string& line)
{
  string auxStr="";
  CFint StchCoef=1,iSpc;
  bool flagProducts=false,readStch=true;
  char currentChar;

  getPartners=false;
  StchVec=0.0;
  idx=0;

  while((currentChar=line[idx])!=':')
  {
    if (currentChar!=' ')
    {
        if (currentChar=='+')
        {
            if(strcoll(auxStr.c_str(),"M")==0) getPartners=true;
            else
            {
              iSpc=0;
              while (strcoll(auxStr.c_str(),_speciesNames[iSpc].c_str())!=0 && iSpc<_NS) iSpc++;
              if (iSpc==_NS)
              {
                CFout << "Wrong format of chemistry file.\n";
                CFout << "At segment: " << auxStr.c_str() <<  "\n";
                abort();
              }

              StchVec[iSpc]=(flagProducts)? StchCoef:-1*StchCoef;
            }

            StchCoef=1;
            auxStr.clear();
            readStch=true;
        }
        else if (currentChar=='=')
        {
            if (flagProducts)
            {
              CFout << "Wrong format of chemistry file.\n";
              CFout << "At segment: " << auxStr.c_str() <<  "\n";
              abort();
            }

            if(strcoll(auxStr.c_str(),"M")==0) getPartners=true;
            else
            {
              iSpc=0;
              while (strcoll(auxStr.c_str(),_speciesNames[iSpc].c_str())!=0 && iSpc<_NS) iSpc++;
              if (iSpc==_NS)
              {
                CFout << "Wrong format of chemistry file.\n";
                CFout << "At segment: " << auxStr.c_str() <<  "\n";
                abort();
              }

              StchVec[iSpc]=-1*StchCoef;
            }

            StchCoef=1;
            auxStr.clear();
            readStch=true;

            flagProducts=true;
        }
        else
        {
            if(readStch&&!isdigit(currentChar))
            {
              if (strcoll(auxStr.c_str(),"")!=0) StchCoef=atoi(auxStr.c_str());

              auxStr.clear();
              readStch=false;
            }
            auxStr.push_back(currentChar);
        }
    }
    idx++;
    if (idx>=(CFint) line.size())
    {
      CFout << "Wrong format of chemistry file.\n";
      CFout << "At segment: " << auxStr.c_str() <<  "\n";
      abort();
    }
  }
  if(strcoll(auxStr.c_str(),"M")!=0)
  {
    iSpc=0;
    while (strcoll(auxStr.c_str(),_speciesNames[iSpc].c_str())!=0 && iSpc<_NS) iSpc++;
    if (iSpc==_NS)
    {
      CFout << "Wrong format of chemistry file.\n";
      CFout << auxStr.c_str() <<  "\n";
      abort();
    }
    StchVec[iSpc]=StchCoef;
  }

  auxStr.resize(0);
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::readChemCoefs(CFint& idx,RealVector& chemCoefs,string& line)
{
    stringstream StrStream;
    string auxStr="";
    CFuint iCoef=0;

    idx++;

    while ((CFuint) idx<line.size())
    {
      if (line[idx] == '/')
      {
        StrStream.str(auxStr);
        if (iCoef>=chemCoefs.size())
        {
            CFout << "Wrong format of chemistry file.\n";
            CFout << line.substr(idx+1) <<  "\n";
            abort();
        }
        if ((StrStream>>chemCoefs[iCoef]).fail())
        {
            CFout << "Wrong format of chemistry file.\n";
            CFout << "At segment: " << auxStr.c_str() << "\n";
            abort();
        }

        StrStream.clear();
        auxStr.clear();
        iCoef++;
      }

      else auxStr.push_back(line[idx]);

      idx++;
    }

    StrStream.str(auxStr);
    if (iCoef!=chemCoefs.size()-1)
    {
        CFout << "Wrong format of chemistry file.\n";
        CFout << line <<  "\n";
        abort();
    }
    if ((StrStream>>chemCoefs[iCoef]).fail())
    {
        CFout << "Wrong format of chemistry file.\n";
        CFout << "At segment: " << auxStr.c_str() << "\n";
        abort();
    }

    StrStream.clear();
    auxStr.resize(0);
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::readPartners(vector<CFreal>& dissPartners,string& line)
{
    string auxStr="";
    CFint idx=0, iSpc;
    char currentChar;

    dissPartners.resize(0);

    while((currentChar=line[idx])!='M')
    {
        idx++;
        if (idx>=((CFint) line.size())+1 || currentChar!=' ')
        {
          CFout << "Wrong format of chemistry file.\n";
          CFout << "At line: " << line <<  "\n";
          abort();
        }
    }
    idx++;
    while((currentChar=line[idx])!='=')
    {
        idx++;
        if (idx>=((CFint) line.size())+1 || currentChar!=' ')
        {
          CFout << "Wrong format of chemistry file.\n";
          CFout << "At line: " << line <<  "\n";
          abort();
        }
    }
    idx++;

    while(idx<(CFint) line.size())
    {
      currentChar=line[idx];
      if (currentChar!=' ')
      {
          if (currentChar==',')
          {
              CFout << "start identifying species \n";
              iSpc=0;
              while (strcoll(auxStr.c_str(),_speciesNames[iSpc].c_str())!=0 && iSpc<_NS) iSpc++;
              if (iSpc==_NS)
              {
                CFout << "Wrong format of chemistry file.\n";
                CFout << "At segment: " << auxStr.c_str() <<  "\n";
                abort();
              }

              dissPartners.push_back(iSpc);

              auxStr.clear();
          }
          else
          {
              auxStr.push_back(currentChar);
          }
      }
      idx++;
    }

    iSpc=0;
    while (strcoll(auxStr.c_str(),_speciesNames[iSpc].c_str())!=0 && iSpc<_NS) iSpc++;
    if (iSpc==_NS)
    {
      CFout << "Wrong format of chemistry file.\n";
      CFout << "At segment: " << auxStr.c_str() <<  "\n";
      abort();
    }

    dissPartners.push_back(iSpc);

    auxStr.resize(0);
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::ComputeOmegav(CFdouble& temperature, RealVector& tVec, CFdouble& pressure, CFdouble& rho, const RealVector& ys, RealVector& omega, RealVector& omegav)
{
    CFdouble tau=0.0;
    CFint i,j,imol;
    omegav[0]=0.0;

    switch(_RelaxID)
    {
    case 0:
      for(i=0;i<_NS;i++)
      {
          _Ss[i]=(_Tets[i]!=-1)? ((_TetsID==0)? 3.5*exp(-_Tets[i]/_Tshk):5000):0;
      }

      for(imol=0;imol<(CFint) _molecIDs.size();imol++)
      {
          i=_molecIDs[imol];
          for(j=0;j<_NS;j++)
          {
              tau+=_massTtmS[j]/exp(_Asr(i,j)*(pow(temperature,-1.0/3)-0.015*pow(_NIUsr(i,j),0.25))-18.42);
          }
          tau*=(pressure*9.6892e-6);
          tau=_massTtmTt/tau;

          tau+=1/(sqrt(8*_Rgas*temperature/(_PI*_mmasses[i]))*1e-21*pow(50000/temperature,2)*rho*_massTtmS[i]*_NA);

          omegav[0]+=ys[i]*(_evibAsterTemp[i]-_evibTemp[i])*pow(abs((_Tshk-tVec[0])/(_Tshk-_Tvsshk)),_Ss[i]-1)/tau;

          tau=0.0;
      }

      omegav[0]*=rho;
      break;

    case 1:
      for(imol=0;imol<(CFint) _molecIDs.size();imol++)
      {
          i=_molecIDs[imol];
          for(j=0;j<_NS;j++)
          {
              tau+=_massTtmS[j]/exp(_Asr(i,j)*(pow(temperature,-1.0/3)-0.015*pow(_NIUsr(i,j),0.25))-18.42);
          }
          tau*=(pressure*9.6892e-6);
          tau=_massTtmTt/tau;
          omegav[0]+=ys[i]*(_evibAsterTemp[i]-_evibTemp[i])/tau;

          tau=0.0;
      }

      omegav[0]*=rho;
      break;

    case 2:
      for(imol=0;imol<(CFint) _molecIDs.size();imol++)
      {
          i=_molecIDs[imol];
          for(j=0;j<_NS;j++)
          {
              tau+=_massTtmS[j]/exp(_Asr(i,j)*(pow(temperature,-1.0/3)-0.015*pow(_NIUsr(i,j),0.25))-18.42);
          }
          tau*=(pressure*9.6892e-6);
          tau=_massTtmTt/tau;

          tau+=1/(sqrt(8*_Rgas*temperature/(_PI*_mmasses[i]))*1e-21*pow(50000/temperature,2)*rho*_massTtmS[i]*_NA);

          omegav[0]+=ys[i]*(_evibAsterTemp[i]-_evibTemp[i])/tau;

          tau=0.0;
      }

      omegav[0]*=rho;
      break;

    case 3:
      for(imol=0;imol<(CFint) _molecIDs.size();imol++)
      {
          i=_molecIDs[imol];
          for(j=0;j<_NS;j++)
          {
              tau+=_massTtmS[j]/(exp(_Asr(i,j)*(pow(temperature,-1.0/3)-0.015*pow(_NIUsr(i,j),0.25))-18.42)*101325+1/(sqrt(8*_NA/(_PI*_NIUsr(i,j)*1e-3*_KB*temperature))*1e-21*pow(50000/temperature,2)));
          }
          tau*=pressure;
          tau=_massTtmTt/tau;

          omegav[0]+=ys[i]*(_evibAsterTemp[i]-_evibTemp[i])/tau;

          tau=0.0;
      }

      omegav[0]*=rho;
      break;
    }

    for(imol=0;imol<(CFint) _molecIDs.size();imol++)
    {
        i=_molecIDs[imol];
        omegav[0]+=_factorOmega*omega[i]*_evibTemp[i];
    }
}

//////////////////////////////////////////////////////////////////////////////

ATDModelLibrary::ChemReact::ChemReact(RealVector& stoichVec, RealVector& chemCoefs,CFuint& chemID,CFuint& TempID)
{
  _stoichVec.resize(stoichVec.size());
  _stoichVec=stoichVec;
  _chemCoefs.resize(chemCoefs.size());
  _chemCoefs=chemCoefs;
  _flagDiss=false;
  _dissPartners.resize(0);
  _chemID=chemID;
  _TempID=TempID;
}

//////////////////////////////////////////////////////////////////////////////

ATDModelLibrary::ChemReact::ChemReact(RealVector& stoichVec, RealVector& chemCoefs, vector<CFreal>& dissPartners,CFuint& chemID,CFuint& TempID)
{
  _stoichVec.resize(stoichVec.size());
  _stoichVec=stoichVec;
  _chemCoefs.resize(chemCoefs.size());
  _chemCoefs=chemCoefs;
  _flagDiss=true;
  _dissPartners.reserve(dissPartners.size());
  for(CFuint i=0;i<dissPartners.size();i++)
  {
    _dissPartners.push_back(dissPartners[i]);
  }
  _chemID=chemID;
  _TempID=TempID;
}

//////////////////////////////////////////////////////////////////////////////

ATDModelLibrary::ChemReact::~ChemReact()
{
    _stoichVec.resize(0);
    _chemCoefs.resize(0);
    _dissPartners.resize(0);
}

//////////////////////////////////////////////////////////////////////////////

ATDModelLibrary::ChemReact::ChemReact()
{
    _stoichVec.resize(0);
    _chemCoefs.resize(0);
    _dissPartners.resize(0);
    _flagDiss=false;
}

//////////////////////////////////////////////////////////////////////////////

void ATDModelLibrary::ChemReact::omegaContribution(CFdouble& temperature, RealVector& tVec,CFdouble& pressure,CFdouble& rho,const RealVector& ys,const RealVector& mmasses,RealVector& omega)
{
    CFdouble fRate=1.0,bRate=-1.0;
    CFuint i;

    for(i=0;i<ys.size();i++)
    {
        if (_stoichVec[i]<0) fRate*=pow(rho*ys[i]/mmasses[i],-_stoichVec[i]);
        else bRate*=pow(rho*ys[i]/mmasses[i],_stoichVec[i]);
    }

    if(_flagDiss)
    {
        switch(_TempID)
        {
        case 0:
            fRate*=Kf(sqrt(temperature*tVec[0]));
            break;
        case 1:
            fRate*=Kf(pow(temperature,0.7)*pow(tVec[0],0.3));
            break;
        }
    }
    else fRate*=Kf(temperature);
    bRate*=Kb(temperature);

    fRate+=bRate;

    if(_flagDiss)
    {
        CFdouble Rate=0.0;
        for(i=0;i<_dissPartners.size();i++)
        {
            Rate+=fRate*rho*ys[_dissPartners[i]]/mmasses[_dissPartners[i]];
        }

        omega+=_stoichVec*Rate;
    }
    else omega+=_stoichVec*fRate;
}

//////////////////////////////////////////////////////////////////////////////

CFdouble ATDModelLibrary::ChemReact::Kf(CFdouble temperature)
{
    return _chemCoefs[0]*pow(temperature,_chemCoefs[1])*exp(-_chemCoefs[2]/temperature);
}

//////////////////////////////////////////////////////////////////////////////

CFdouble ATDModelLibrary::ChemReact::Kb(CFdouble temperature)
{
    CFdouble Kback=0.0,z=10000/temperature;
    CFuint i;

    switch(_chemID)
    {
    case 0:
        Kback=_chemCoefs[7];
        for(i=6;i>2;i--)
        {
            Kback*=z;
            Kback+=_chemCoefs[i];
        }

        Kback=(Kback<-100)? -100:Kback;

        Kback=Kf(temperature)*exp(-Kback);
        break;
    case 1:
        for(i=7;i>4;i--)
        {
            Kback+=_chemCoefs[i];
            Kback*=z;
        }
        Kback+=_chemCoefs[4]*log(z);
        Kback+=_chemCoefs[3];

        Kback=(Kback<-100)? -100:Kback;

        Kback=Kf(temperature)*exp(-Kback);
        break;
    case 2:
        for(i=7;i>5;i--)
        {
            Kback+=_chemCoefs[i];
            Kback*=z;
        }
        Kback+=_chemCoefs[5]*log(z);
        Kback+=_chemCoefs[4];
        Kback+=_chemCoefs[3]/z;
        Kback=(Kback<-100)? -100:Kback;

        Kback=Kf(temperature)*exp(-Kback);
        break;
    case 3:
        Kback=_chemCoefs[3]*pow(temperature,_chemCoefs[4])*exp(-_chemCoefs[5]/temperature);
        break;
    }

    return Kback;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace ATDModel

} // namespace Physics

} // namespace COOLFluiDs

//////////////////////////////////////////////////////////////////////////////

