#include "ATDModel/CombustionModelLibrary.hh"
#include "ATDModel/ATDModel.hh"
#include "Common/CFLog.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/BadValueException.hh"
#include "Framework/PhysicalModel.hh"
#include "Common/Stopwatch.hh"
#include "Environment/DirPaths.hh"
#include "Common/OSystem.hh"
#include "Common/PE.hh"
#include "Common/StringOps.hh"
#include "CombustionModelLibrary.hh"

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

Environment::ObjectProvider<CombustionModelLibrary,
                            PhysicalPropertyLibrary,
                            ATDModelModule,
                            1>
combustionModelLibraryProvider("CombustionModel");

//////////////////////////////////////////////////////////////////////////////

void CombustionModelLibrary::defineConfigOptions(Config::OptionList& options)
{
}

//////////////////////////////////////////////////////////////////////////////

CombustionModelLibrary::CombustionModelLibrary(const std::string& name)
  : ATDModelLibrary(name),
    _xs(), 
    _AsMu1(), _BsMu1(), _CsMu1(), _DsMu1(), _AsMu2(), _BsMu2(), _CsMu2(), _DsMu2(), 
    _AsK1(), _BsK1(), _CsK1(), _DsK1(), _AsK2(), _BsK2(), _CsK2(), _DsK2(), 
    _a1cp1(), _a2cp1(), _a3cp1(), _a4cp1(), _a5cp1(), 
    _a1cp2(),_a2cp2(),_a3cp2(),_a4cp2(),_a5cp2(),
    _a1h1(),_a2h1(),_a3h1(),_a4h1(),_a5h1(),_a6h1(),
    _a1h2(),_a2h2(),_a3h2(),_a4h2(),_a5h2(),_a6h2(),
    _KappasTemp(),
    _CpTemp(),
    _hsTemp(),
    _hsP(),
    _omegaTempDot()
{
  addConfigOptionsTo(this);
  
  _mixtureName = "";
  _thermoName = "polynomial";
  _chemName = "";
}

//////////////////////////////////////////////////////////////////////////////

CombustionModelLibrary::~CombustionModelLibrary()
{
}

//////////////////////////////////////////////////////////////////////////////

void CombustionModelLibrary::configure ( Config::ConfigArgs& args )
{
  ATDModelLibrary::configure(args);

  // here you are reading in the input file settings
}

//////////////////////////////////////////////////////////////////////////////

void CombustionModelLibrary::unsetup()
{ 
  if(isSetup()) {
    
    // De-allocate local working arrays
    _xs.resize(0);
    _AsMu1.resize(0);
    _BsMu1.resize(0);
    _CsMu1.resize(0);
    _DsMu1.resize(0);
    _AsMu2.resize(0);
    _BsMu2.resize(0);
    _CsMu2.resize(0);
    _DsMu2.resize(0);
    _AsK1.resize(0);
    _BsK1.resize(0);
    _CsK1.resize(0);
    _DsK1.resize(0);
    _AsK2.resize(0);
    _BsK2.resize(0);
    _CsK2.resize(0);
    _DsK2.resize(0);
    _a1cp1.resize(0);
    _a2cp1.resize(0);
    _a3cp1.resize(0);
    _a4cp1.resize(0);
    _a5cp1.resize(0);
    _a1cp2.resize(0);
    _a2cp2.resize(0);
    _a3cp2.resize(0);
    _a4cp2.resize(0);
    _a5cp2.resize(0);
    _a1h1.resize(0);
    _a2h1.resize(0);
    _a3h1.resize(0);
    _a4h1.resize(0);
    _a5h1.resize(0);
    _a6h1.resize(0);
    _a1h2.resize(0);
    _a2h2.resize(0);
    _a3h2.resize(0);
    _a4h2.resize(0);
    _a5h2.resize(0);
    _a6h2.resize(0);
    _KappasTemp.resize(0);
    _CpTemp.resize(0);
    _hsTemp.resize(0);
    _hsP.resize(0);
    _omegaTempDot.resize(0);
    
    ATDModelLibrary::unsetup();
  }
}

//////////////////////////////////////////////////////////////////////////////

void CombustionModelLibrary::setMoleculesIDs(std::vector<CFuint>& v)
{
//  // is this generic enough ???
////   CFout << "Start of setMoleculesIDs"<< "\n";
//  v.reserve(_NS);
//  for (CFint i = 0; i < _NS; ++i) {
//    if (_CA[i] > 1) {
//      v.push_back(i);
//    }
//  }
}

//////////////////////////////////////////////////////////////////////////////

CFdouble CombustionModelLibrary::eta(CFdouble& temp, CFdouble& pressure, CFreal* tVec)
{ 
  //CFout <<"Entering eta() function \n";
  
  //CFout << "Start of eta function" << "\n";
  CFdouble eta = 0.0, PHIs;
 
  
 
  
  _EpsKb = 337.6241;
  _Tstar = 0.003*temp;
  _Ov = 1.16145*std::pow(_Tstar, -0.14874)+0.52487*exp(-0.77320*_Tstar)+2.16178*exp(-2.43787*_Tstar);
  _Fc = 1-0.2756*_AF;
  
  // C4H6 Butadiene requires specific viscosity model; it is placed as first member of species viscosity
  // array (place [0]) - _mmasses[0] is converted inside the calculation from kg to g again
  _NIUsTemp[0]=40.785*((_Fc*sqrt(_mmasses[0]*1000*temp))/(std::pow(_Vc,0.6667)*_Ov));
  _NIUsTemp[0]=_NIUsTemp[0]*std::pow(10.,-7.);
   // CFout << _NIUsTemp[0] << "\n";
  
  //CFout << "Eta function temperature " << temp << "\n";
  //CFout << "Eta function _NIUsTemp[0] " << _NIUsTemp[0] << "\n";
  

  _massTtmS[0]=_ys[0]/_mmasses[0];
  
  
  switch(_chemID) 
  {
  //CFout << "Temperature in eta" << temp << "\n";  
  case 0:
  //CFout <<"Switching to eta() function case 0 \n";
  CFint i;//,j;
  for(i=1;i< _NS;i++)
    {
        _massTtmS[i]=_ys[i]/_mmasses[i];
	if (i==1){
	  if (temp <= 373.2) {
	     _NIUsTemp[i]=1.2277*std::pow(10.,-5.);
// 	     ofstream viscfile((_libPath+"validation/visc1.txt").c_str());
// 	     viscfile.open("visc1.txt");
// 	     CFout << _NIUsTemp[i] << "\n";
// 	     viscfile << "prova";
// 	     viscfile.close();
	     
	  }
	  else if (temp>373.2 && temp<1075.0) {
	          _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]);
		  _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
// 		   CFout << _NIUsTemp[i] << "\n";
	  }
	  else {
	       _NIUsTemp[i]=exp(_AsMu2[i]*log(temp)+_BsMu2[i]/temp+_CsMu2[i]/(temp*temp)+_DsMu2[i]);
	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
// 	        CFout << _NIUsTemp[i] << "\n";
	  } 
	}
	else if (i != 1) {
	   if (temp < 1000) {
	      _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]);
	      _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
// 	       CFout << _NIUsTemp[i] << "\n";
	  }
	  else {
	       _NIUsTemp[i]=exp(_AsMu2[i]*log(temp)+_BsMu2[i]/temp+_CsMu2[i]/(temp*temp)+_DsMu2[i]);
	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
// 	        CFout << _NIUsTemp[i] << "\n";
	  }
	}
      }

  
  break;
      
  case 1:
  //CFout << "Temperature in eta" << temp << "\n";  
  
  //CFout <<"Switching to eta() function case 1 \n";
 
  for(i=1;i< _NS;i++)
    {
        _massTtmS[i]=_ys[i]/_mmasses[i];
	//H2O HAS to be the species in position [1] of the vectors (.cmix file accordingly!)
	if (i==1){
	  if (temp <= 373.2) {
	     _NIUsTemp[i]=1.2277*std::pow(10.,-5.);
	  }
	    else if (temp>373.2 && temp<=1075.0) {
	    _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]);
	    _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	  else {
	       _NIUsTemp[i]=exp(_AsMu2[i]*log(temp)+_BsMu2[i]/temp+_CsMu2[i]/(temp*temp)+_DsMu2[i]);
	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  } 
	}
	else if (i != 1) {
	   if (temp < 1000) {
	      _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]);
	      _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	  else {
	       _NIUsTemp[i]=exp(_AsMu2[i]*log(temp)+_BsMu2[i]/temp+_CsMu2[i]/(temp*temp)+_DsMu2[i]);
	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	}
      }
  // C4H6 Butadiene requires specific viscosity model; it is placed as first member of species viscosity
  // array (place [0])

//   CFout <<"_NIUsTemp from eta"<< _NIUsTemp << "\n";
  
  break;
  
  case 2:
    
    //CFout << "Temperature in eta" << temp << "\n";
    
  //CFout <<"Switching to eta() function case 2 \n";
  for(i=1;i< _NS;i++)
    {
        _massTtmS[i]=_ys[i]/_mmasses[i];
	//H2O HAS to be the species in position [1] of the vectors (.cmix file accordingly!)
	if (i==1){
	  if (temp <= 373.2) {
	     _NIUsTemp[i]=1.2277*std::pow(10.,-5.);
	  }
	    else if (temp>373.2 && temp<=1075.0) {
	    _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]);
	    _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	  else {
	       _NIUsTemp[i]=exp(_AsMu2[i]*log(temp)+_BsMu2[i]/temp+_CsMu2[i]/(temp*temp)+_DsMu2[i]);
	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  } 
	}
	//H HAS to be the species in position [2] of the vectors (.cmix file accordingly!)
	else if (i==2) {
	  if (temp <=1000.0) {
	     _NIUsTemp[i]=1.4236*std::pow(10.,-5.);
	  }  
	  else {
	       _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]); 
	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	}
	  //O HAS to be the species in position [2] of the vectors (.cmix file accordingly!)
	else if (i==3) {
	  if (temp <=1000.0) {
	     _NIUsTemp[i]=4.9965*std::pow(10.,-5.);
	  }
	  else {
	       _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]); 
	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	}
	  // This can probably be substituted by a simple "else" without if
	else { //if (i != 1 && i != 2 && i != 3) {
	   if (temp < 1000) {
	      _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]); 
	      _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	
	  else {
	       _NIUsTemp[i]=exp(_AsMu2[i]*log(temp)+_BsMu2[i]/temp+_CsMu2[i]/(temp*temp)+_DsMu2[i]);
	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	}
      }
  // C4H6 Butadiene requires specific viscosity model; it is placed as first member of species viscosity
  // array (place [0])

   break;
    
   }
   
 
CFint i,j;
    // Mixture viscosity calculation as sum weighted over PHI.
    for(i=0;i< _NS;i++)
    {
        PHIs=0.0;
	// PHIs is the same for ATDModelLibrary and CombustionModelLibrary
        for(j=0;j<_NS;j++)
        {
            PHIs+=_massTtmS[j]*std::pow(1+sqrt(_NIUsTemp[i]/_NIUsTemp[j])*std::pow(_mmasses[j]/_mmasses[i],0.25),2)/sqrt(8*(1+_mmasses[i]/_mmasses[j]));
        }

        eta+=_massTtmS[i]*_NIUsTemp[i]/PHIs;
    }

  //CFout << "End of eta function" << "eta: " << eta << "\n";
  return eta;
  
}

//////////////////////////////////////////////////////////////////////////////

CFdouble CombustionModelLibrary::lambdaNEQ(CFdouble& temp, CFdouble& pressure)//,
                                           //CFreal* tVec)
// Before it was (CFdouble& temperatureerature,CFdouble& pressure)
{ 
 // CFout << "Start of lambdaNEQ function" << "\n";
  
  CFdouble lambda = 0.0, PHIs;
  
  // Alpha factor for C4H6 thermal conductivity
  CFdouble _alphaK = 0.0;
  // Beta factor for C4H6 thermal conductivity
  CFdouble _betaK = 0.0;
  // Zeta factor for C4H6 thermal conductivity
  CFdouble _zetaK = 0.0;
  // Psi factor for C4H6 thermal conductivity
  CFdouble _PSIK = 0.0;
  
  //   Function lambdaNEQ has to calculate thermal conductivity of each species given temperatureerature
  //   viscosiy-like method.
  //   In order to do this, viscosity of each species and of mixture has to be re-calculate in order to calculate
  //   weight function PHIs
  //   Then the same method has to be applied changnin NIUTemp to LambdaTemp with AsK, BsK, CsK, DsK
  _EpsKb = 337.6241;
  _Tstar = 0.003*temp;
  _Ov = 1.16145*std::pow(_Tstar, -0.14874)+0.52487*exp(-0.77320*_Tstar)+2.16178*exp(-2.43787*_Tstar);
  _Fc = 1-0.2756*_AF;
  
  // C4H6 Butadiene requires specific viscosity model; it is placed as first member of species viscosity
  // array (place [0]) - _mmasses[0] is converted inside the calculation from kg to g again
  _NIUsTemp[0]=40.785*((_Fc*sqrt(_mmasses[0]*1000*temp))/(std::pow(_Vc,0.6667)*_Ov));
  _NIUsTemp[0]=_NIUsTemp[0]*std::pow(10.,-7.);
  _massTtmS[0]=_ys[0]/_mmasses[0];
 
// CFout << "temp" << temp << "\n";
//  CFout << "temperature" << temperature << "\n";
//  CFout << "_mmasses[0]" << _mmasses[0] << "\n";
//  CFout << "NIUsTemp from lambdaNEQ: " << _NIUsTemp[0] << "\n";
  

  
// Viscosity is necessary for thermal conductivity weighing calculation (through PHIs)
// and it is here re-calculated
  switch(_chemID)
  {
    
    //CFout << "Temperature in lambdaNEQ" << temp << "\n";
  case 0:
  
  CFint i;//,j;
  for(i=1;i< _NS;i++)
    {
        _massTtmS[i]=_ys[i]/_mmasses[i];
	if (i==1){
	  if (temp <= 373.2) {
	     _NIUsTemp[i]=1.2277*std::pow(10.,-5.);
	  }
	  else if (temp>373.2 && temp<1075.0) {
	          _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]);
		  _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	  else {
	       _NIUsTemp[i]=exp(_AsMu2[i]*log(temp)+_BsMu2[i]/temp+_CsMu2[i]/(temp*temp)+_DsMu2[i]);
	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  } 
	}
	else if (i != 1) {
	   if (temp < 1000) {
	      _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]); 
	      _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	  else {
	       _NIUsTemp[i]=exp(_AsMu2[i]*log(temp)+_BsMu2[i]/temp+_CsMu2[i]/(temp*temp)+_DsMu2[i]); 
	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	}
      }
//    CFout << "NIUsTemp from lambdaNEQ: " << _NIUsTemp << "\n";  
  break;
      
  case 1:
  
 // CFout << "Temperature in lambdaNEQ" << temp << "\n";
  for(i=1;i< _NS;i++)
    {
        _massTtmS[i]=_ys[i]/_mmasses[i];
	//H2O HAS to be the species in position [1] of the vectors (.cmix file accordingly!)
	if (i==1){
	  if (temp <= 373.2) {
	     _NIUsTemp[i]=1.2277*std::pow(10.,-5.);
	  }
	    else if (temp>373.2 && temp<=1075.0) {
	    _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]);
	    _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	  else {
	       _NIUsTemp[i]=exp(_AsMu2[i]*log(temp)+_BsMu2[i]/temp+_CsMu2[i]/(temp*temp)+_DsMu2[i]);
	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  } 
	}
	else if (i != 1 ) {
	   if (temp < 1000) {
	      _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]); 
	      _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	  else {
	       _NIUsTemp[i]=exp(_AsMu2[i]*log(temp)+_BsMu2[i]/temp+_CsMu2[i]/(temp*temp)+_DsMu2[i]); 
	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	}
      }

//   CFout << "NIUsTemp from lambdaNEQ: " << _NIUsTemp << "\n";
  
  break;
  
  case 2:
  
 // CFout << "Temperature in lambdaNEQ" << temp << "\n";  

  for(i=1;i< _NS;i++)
    {
        _massTtmS[i]=_ys[i]/_mmasses[i];
	//H2O HAS to be the species in position [1] of the vectors (.cmix file accordingly!)
	if (i==1){
	  if (temp <= 373.2) {
	     _NIUsTemp[i]=1.2277*std::pow(10.,-5.);
	  }
	    else if (temp>373.2 && temp<=1075.0) {
	    _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]);
	    _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	  else {
	       _NIUsTemp[i]=exp(_AsMu2[i]*log(temp)+_BsMu2[i]/temp+_CsMu2[i]/(temp*temp)+_DsMu2[i]);
	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  } 
	}
	//H HAS to be the species in position [2] of the vectors (.cmix file accordingly!)
	else if (i==2) {
	  if (temp <=1000.0) {
	     _NIUsTemp[i]=1.4236*std::pow(10.,-5.);
	  }  
	  else {
	       _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]); 
	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	}
	//O HAS to be the species in position [2] of the vectors (.cmix file accordingly!)
	else if (i==3) {
	  if (temp <=1000.0) {
	     _NIUsTemp[i]=4.9965*std::pow(10.,-5.);
	  }  
	  else {
	       _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]); 
	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	}
	// This can probably be substituted by a simple "else" without if
	else if (i != 1 && i != 2 && i != 3) {
	   if (temp < 1000) {
	      _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]);
	      _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	  else {
	       _NIUsTemp[i]=exp(_AsMu2[i]*log(temp)+_BsMu2[i]/temp+_CsMu2[i]/(temp*temp)+_DsMu2[i]); 
	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	}
      }
  
   break;
    
   }
    
 
// CALCULATION OF THERMAL CONDUCTIVITY 
  // C4H6 Butadiene requires specific thermal conductivity model; it is placed as first member of species thermal
  // array (place [0]). It requires also computation of specific heat as first step.

if (temp <= 1000) {
	       _CpTemp[0] = _Rgas * (_a1cp1[0] + (_a2cp1[0] * temp) + (_a3cp1[0] * std::pow(temp,2)) + (_a4cp1[0] * std::pow(temp,3)) + (_a5cp1[0] * std::pow(temp,4))); 
		 // old version for CEAS-Thermobuild coeffs, in data_old folder
	     //_CpTemp[0]= _a1cp1[0]*std::pow(temp,4)+_a2cp1[0]*std::pow(temp,3)+_a3cp1[0]*std::pow(temp,2)+_a4cp1[0]*temp+_a5cp1[0];
	  }
	  else {
	       _CpTemp[0] = _Rgas * (_a1cp2[0] + (_a2cp2[0] * temp) + (_a3cp2[0] * std::pow(temp,2)) + (_a4cp2[0] * std::pow(temp,3)) + (_a5cp2[0] * std::pow(temp,4))); 
		  // old version for CEAS-Thermobuild coeffs, in data_old folder 
		  //_CpTemp[0]= _a1cp2[0]*std::pow(temp,4)+_a2cp2[0]*std::pow(temp,3)+_a3cp2[0]*std::pow(temp,2)+_a4cp2[0]*temp+_a5cp2[0];
	  } 
      
  //_CpC4H6 = _CpTemp[0];
  //CFout << "_Rgas:" << _Rgas << "\n";
  //CFout << _CpTemp[0] << "\n";
  
  _AF = 0.192;
  
  _alphaK = ((_CpTemp[0]-_Rgas)/_Rgas)-1.5;
  //CFout << "_alphaK:" << _alphaK << "\n";
  //CFout << _CpTemp[0];
  _betaK = 0.7862 - 0.7109*_AF + 1.3168 * (_AF*_AF);
  //CFout << "_betaK:" << _betaK << "\n";
  _Tr = temp/_Tc;
  _zetaK = 2.0 + 10.5*(_Tr*_Tr);
  //CFout << "_zetaK:" << _zetaK << "\n";
  _PSIK = 1 + _alphaK*((0.215+0.28288*_alphaK-1.061*_betaK+0.26665*_zetaK)/(0.6366+_betaK*_zetaK+1.061*_alphaK*_betaK));
  //CFout << "_PSIK:" << _PSIK << "\n";
  
  
  _KappasTemp[0]= ((3.75*_PSIK)/((_CpTemp[0]/_Rgas)-1))*((_NIUsTemp[0]*(_CpTemp[0]-_Rgas))/(_mmasses[0]));
  
  _massTtmS[0]=_ys[0]/_mmasses[0];
  

 switch(_chemID)
  {
  case 0:
  
  CFint i;
  for(i=1;i< _NS;i++)
    {
        _massTtmS[i]=_ys[i]/_mmasses[i];
	if (i==1){
	  if (temp <= 373.2) {
	     _KappasTemp[i]=0.02499;
	     
	  }
	  else if (temp>373.2 && temp<1075.0) {
	          _KappasTemp[i]=exp(_AsK1[i]*log(temp)+_BsK1[i]/temp+_CsK1[i]/(temp*temp)+_DsK1[i]);
		  _KappasTemp[i]=_KappasTemp[i]*std::pow(10.,-4);
	  }
	  else {
	       _KappasTemp[i]=exp(_AsK2[i]*log(temp)+_BsK2[i]/temp+_CsK2[i]/(temp*temp)+_DsK2[i]);
	       _KappasTemp[i]=_KappasTemp[i]*std::pow(10.,-4);
	  } 
	}
	else if (i != 1) {
	   if (temp < 1000) {
	      _KappasTemp[i]=exp(_AsK1[i]*log(temp)+_BsK1[i]/temp+_CsK1[i]/(temp*temp)+_DsK1[i]);
	      _KappasTemp[i]=_KappasTemp[i]*std::pow(10.,-4);
	  }
	  else {
	       _KappasTemp[i]=exp(_AsK2[i]*log(temp)+_BsK2[i]/temp+_CsK2[i]/(temp*temp)+_DsK2[i]);
	       _KappasTemp[i]=_KappasTemp[i]*std::pow(10.,-4);
	  }
	}
      }
      

  
  break;
      
  case 1:
    

  for(i=1;i< _NS;i++)
    {
        _massTtmS[i]=_ys[i]/_mmasses[i];
	//H2O HAS to be the species in position [1] of the vectors (.cmix file accordingly!)
	if (i==1){
	  if (temp <= 373.2) {
	     _KappasTemp[i]=0.02499;
	  }
	    else if (temp>373.2 && temp<=1075.0) {
	    _KappasTemp[i]=exp(_AsK1[i]*log(temp)+_BsK1[i]/temp+_CsK1[i]/(temp*temp)+_DsK1[i]);
	    _KappasTemp[i]=_KappasTemp[i]*std::pow(10.,-4);
	  }
	  else {
	       _KappasTemp[i]=exp(_AsK2[i]*log(temp)+_BsK2[i]/temp+_CsK2[i]/(temp*temp)+_DsK2[i]);
	       _KappasTemp[i]=_KappasTemp[i]*std::pow(10.,-4);
	  } 
	}
	else if (i != 1) {
	   if (temp < 1000) {
	      _KappasTemp[i]=exp(_AsK1[i]*log(temp)+_BsK1[i]/temp+_CsK1[i]/(temp*temp)+_DsK1[i]);
	      _KappasTemp[i]=_KappasTemp[i]*std::pow(10.,-4);
	  }
	  else {
	       _KappasTemp[i]=exp(_AsK2[i]*log(temp)+_BsK2[i]/temp+_CsK2[i]/(temp*temp)+_DsK2[i]);
	       _KappasTemp[i]=_KappasTemp[i]*std::pow(10.,-4);
	  }
	}
      }
//   CFout <<"KappasTemp: " << _KappasTemp << "\n";
  
  break;
  
  case 2:
   

  for(i=1;i< _NS;i++)
    {
        _massTtmS[i]=_ys[i]/_mmasses[i];
	//H2O HAS to be the species in position [1] of the vectors (.cmix file accordingly!)
	if (i==1){
	  if (temp <= 373.2) {
	     _KappasTemp[i]=0.02499;
	  }
	    else if (temp>373.2 && temp<=1075.0) {
	    _KappasTemp[i]=exp(_AsK1[i]*log(temp)+_BsK1[i]/temp+_CsK1[i]/(temp*temp)+_DsK1[i]);
	    _KappasTemp[i]=_KappasTemp[i]*std::pow(10.,-4);
	  }
	  else {
	       _KappasTemp[i]=exp(_AsK2[i]*log(temp)+_BsK2[i]/temp+_CsK2[i]/(temp*temp)+_DsK2[i]);
	       _KappasTemp[i]=_KappasTemp[i]*std::pow(10.,-4);
	  } 
	}
	//H HAS to be the species in position [2] of the vectors (.cmix file accordingly!)
	else if (i==2) {
	  if (temp <=1000.0) {
	     _KappasTemp[i]=0.440387;
	  }  
	  else {
	       _KappasTemp[i]=exp(_AsK1[i]*log(temp)+_BsK1[i]/temp+_CsK1[i]/(temp*temp)+_DsK1[i]); 
	       _KappasTemp[i]=_KappasTemp[i]*std::pow(10.,-4);
	  }
	}
	//O HAS to be the species in position [2] of the vectors (.cmix file accordingly!)
	else if (i==3) {
	  if (temp <=1000.0) {
	     _KappasTemp[i]=0.097372;
	  }  
	  else {
	       _KappasTemp[i]=exp(_AsK1[i]*log(temp)+_BsK1[i]/temp+_CsK1[i]/(temp*temp)+_DsK1[i]);
	       _KappasTemp[i]=_KappasTemp[i]*std::pow(10.,-4);
	  }
	}
	// This can probably be substituted by a simple "else" without if
	else if (i != 1 && i != 2 && i != 3) {
	   if (temp < 1000) {
	      _KappasTemp[i]=exp(_AsK1[i]*log(temp)+_BsK1[i]/temp+_CsK1[i]/(temp*temp)+_DsK1[i]);
	      _KappasTemp[i]=_KappasTemp[i]*std::pow(10.,-4);
	  }
	  else {
	       _KappasTemp[i]=exp(_AsK2[i]*log(temp)+_BsK2[i]/temp+_CsK2[i]/(temp*temp)+_DsK2[i]);
	       _KappasTemp[i]=_KappasTemp[i]*std::pow(10.,-4);
	  }
	}
      }
    
   
   
  
   break;
    
   }
    CFint i,j;
    for(i=0;i< _NS;i++)
    {
        PHIs=0.0;
	// PHIs is the same for ATDModelLibrary and CombustionModelLibrary
        for(j=0;j<_NS;j++)
        {
            PHIs+=_massTtmS[j]*std::pow(1+sqrt(_NIUsTemp[i]/_NIUsTemp[j])*std::pow(_mmasses[j]/_mmasses[i],0.25),2)/sqrt(8*(1+_mmasses[i]/_mmasses[j]));
        }

        lambda+=_massTtmS[i]*_KappasTemp[i]/PHIs;
    }
	
	//CFout << "End of lambdaNEQ function -> lambda " << lambda << "\n";
  return lambda;
 // CFout << "End of lambdaNEQ function" << "\n";

 }
 
// //   CFint i,j;
// //   for(i=1;i< _NS;i++)
// //     {
// //         _massTtmS[i]=_ys[i]/_mmasses[i];
// // 	//added "_Ds"; it has to be added to input files, to be declared and to be read from file
// // 	
// // 	//Insert "if" here to sort between temp<1000 and temp>1000 for all except O, H, H2O, C4H6
// // 	//Insert another "if" accounting for H2O use of 1075 instead of 1000
// // 	//H2O must have constant value from 300 to 373.2
// // 	//O and H must have constant value before 1000
// // 	//_NIUsTemp[i]=exp(_AsMu[i]*log(temp)+_BsMu[i]/temp+_CsMu[i]/(temp*temp)+_DsMu[i]);
// // 	//_LambdaTemp yet to declare at the right place 
// // 	//_LambdaTemp[i]=exp(_AsK[i]*log(temp)+_BsK[i]/temp+_CsK[i]/(temp*temp)+_DsK[i]);
// //     }

//////////////////////////////////////////////////////////////////////////////

CFdouble CombustionModelLibrary::CpPoly(CFdouble& temp)//,CFint& i)
{
    //CFout << "Inside CpPoly function" << "\n";   
   
    CFdouble CpP = 0.0;
    
    
switch(_thermoID)
  {
  case 0: //Polynomial interpolation of Cp

  //CFout << "Inside Thermo-Switch" << "\n";
  
  CFint i;
  for(i=0;i< _NS;i++)
    {   
		if (temp <= 1000) {
			_CpTemp[i] = _Rgas * (_a1cp1[i] + (_a2cp1[i] * temp) + (_a3cp1[i] * std::pow(temp,2)) + (_a4cp1[i] * std::pow(temp,3)) + (_a5cp1[i] * std::pow(temp,4))); 
			
			
		}
		else {
			_CpTemp[i] = _Rgas * (_a1cp2[i] + (_a2cp2[i] * temp) + (_a3cp2[i] * std::pow(temp,2)) + (_a4cp2[i] * std::pow(temp,3)) + (_a5cp2[i] * std::pow(temp,4))); 
			
		}
		
	 // old version for CEAS-Thermobuild coeffs, in data_old folder	
      // 	if (temp <= 1000) {
	  //   _CpTemp[i]= _a1cp1[i]*std::pow(temp,4)+_a2cp1[i]*std::pow(temp,3)+_a3cp1[i]*std::pow(temp,2)+_a4cp1[i]*temp+_a5cp1[i];
	  //}
	  //else {
	  //     _CpTemp[i]= _a1cp2[i]*std::pow(temp,4)+_a2cp2[i]*std::pow(temp,3)+_a3cp2[i]*std::pow(temp,2)+_a4cp2[i]*temp+_a5cp2[i];      
	  //}
	 
	 // Specific heat from J/mol to J/kg
	 _CpTemp[i] /=_mmasses[i];
	 
	 CpP += _CpTemp[i]*_ys[i];
    }
    
//       for(i=0;i< _NS;i++)
//     {   
//       CpP += _CpTemp[i]*_ys[i];
//       
//     }
//     
//  CFdouble massTot=0.0;
// 
//   for (CFint is = 0; is < _NS; ++is) {
//     cf_always_assert (_ys[is] > 1.1);
// 
//     const CFreal mm = _ys[is]/_mmasses[is];
//     massTot += mm;
//     _xs[is] = mm;
//   }
//   
//      _xs *= 1./massTot;
// 	 
// 
//   for(i=0;i< _NS;i++)
//     {   
//       CpP += _CpTemp[i]*_xs[i];
//       
//     }

// Mixture specific heat CpP weighted on mass fraction



 
  //CFout << "Before CpPoly Return." << "\n";
  

  }
return CpP;
}
    
    
    
//     if (temp<=200){
// //        CFout << "Warning: Temperature outside allowable range, " << temp << "K is less than 200K ( at enthalpyTeqOverRMB() ). Extrapolating. \n";
//         iind=i*_numTC;
//     }
//     else if (temp<=1000) iind=i*_numTC;
//     else if (temp<=6000) iind=i*_numTC+_numTCR;
//     else if (temp<=20000) iind=i*_numTC+2*_numTCR;
//     else {
// //        CFout << "Warning: Temperature outside allowable range, " << temp << "K is greater than 20000K ( at enthalpyTeqOverRMB() ). Extrapolating. \n";
//         iind=i*_numTC+2*_numTCR;
//     }
// 
//     for (CFint j = 5; j >= 1; --j) {
//       CpOR *= temp;
//       CpOR += _thermoCoefs[iind+j+1];
//     }
// 
//     CpOR += (_thermoCoefs[iind]/temp+_thermoCoefs[iind+1])/temp;
// 
//     return CpOR/_mmasses[i];
// }

//////////////////////////////////////////////////////////////////////////////

void CombustionModelLibrary::EnthPoly(CFdouble& temp,CFint& i)// CFdouble& _hsP)
{
    //CFint iind;

//CFout << "Start of EnthPoly function" << "\n";    
    
switch(_thermoID)
  {
    
  //CFout << "Temperature in EnthPoly" << temp << "\n";  
  case 0: //Polynomial interpolation of Cp

  CFint i;//,j;
  for(i=0;i< _NS;i++)
    {
		if (temp <= 1000) {
			_hsTemp[i]= _Rgas*temp*((_a1h1[i]) + (_a2h1[i] * temp)/2 + (_a3h1[i] * std::pow(temp,2))/3 + (_a4h1[i] * std::pow(temp,3))/4 + (_a5h1[i] * std::pow(temp,4))/5 + (_a6h1[i]/temp)); 
			// old version for CEAS-Thermobuild coeffs, in data_old folder
			//_hsTemp[i]= _a1h1[i]*std::pow(temp,5)+_a2h1[i]*std::pow(temp,4)+_a3h1[i]*std::pow(temp,3)+_a4h1[i]*std::pow(temp,2)+_a5h1[i]*temp+_a6h1[i];
		}
		else {
			_hsTemp[i]= _Rgas*temp*((_a1h2[i]) + (_a2h2[i] * temp)/2 + (_a3h2[i] * std::pow(temp,2))/3 + (_a4h2[i] * std::pow(temp,3))/4 + (_a5h2[i] * std::pow(temp,4))/5 + (_a6h2[i]/temp));
			// old version for CEAS-Thermobuild coeffs, in data_old folder
			//_hsTemp[i]= _a1h2[i]*std::pow(temp,5)+_a2h2[i]*std::pow(temp,4)+_a3h2[i]*std::pow(temp,3)+_a4h2[i]*std::pow(temp,2)+_a5h2[i]*temp+_a6h2[i];;
		} 
	}
	  // Conversion from "per mole unit" to "per mass unit"
	 _hsTemp[i] /= _mmasses[i];
	  
	  
   //   CFout << "Before EnthPoly last line." << "\n";
         _hsP[i] = _hsTemp[i];
	 
	 
//         _massTtmS[i]=_ys[i]/_mmasses[i];
// 	
// 	CFdouble massTot=0.0;
// 
// 	for (CFint is = 0; is < _NS; ++is) {
// 	cf_always_assert (_ys[is] > 1.1);
// 
// 	const CFreal mm = _ys[is]/_mmasses[is];
// 	massTot += mm;
// 	_xs[is] = mm;
// 	}
// 	
//     _xs *= 1./massTot;
// 	
// 	_meanmass +=_mmasses[i];
// 	_meanmass = _meanmass/_NS;
// 	_xs[i] = _ys[i]*_meanmass;
// 	_xs[i] = _xs[i]/_mmasses[i];
// 	

      
    }
  
  }





//////////////////////////////////////////////////////////////////////////////

void CombustionModelLibrary::frozenGammaAndSoundSpeed(CFdouble& temp,
						      CFdouble& pressure,
						      CFdouble& rho,
						      CFdouble& gamma,
						      CFdouble& soundSpeed,
						      RealVector* tVec)
{
  // 
  CFdouble RoverMtot = 0.0;
  CFdouble CvP = 0.0;
  CFdouble CpP = 0.0;
  _massTtmTt = 0.0;
  
//   _meanmass = 0.0;
  
  //CFout << "temp in frozenGammaAndSoundSpeed" << temp << "\n";
  CFint i;
  //_CvTrOR=0.0;
//CFout << "inside frozenGammaAndSoundSpeed" << "\n";
  for ( i = 0; i < _NS; ++i) {
    
     _massTtmS[i] = _ys[i]/_mmasses[i];
     _massTtmTt += _massTtmS[i];
//     
	  if (temp <= 1000) {
		  _CpTemp[i] = _Rgas * (_a1cp1[i] + (_a2cp1[i] * temp) + (_a3cp1[i] * std::pow(temp,2)) + (_a4cp1[i] * std::pow(temp,3)) + (_a5cp1[i] * std::pow(temp,4))); 
		  
		  
	  }
	  else {
		  _CpTemp[i] = _Rgas * (_a1cp2[i] + (_a2cp2[i] * temp) + (_a3cp2[i] * std::pow(temp,2)) + (_a4cp2[i] * std::pow(temp,3)) + (_a5cp2[i] * std::pow(temp,4))); 
		  
	  }
	  
	  // old version for CEAS-Thermobuild coeffs, in data_old folder	
      // 	if (temp <= 1000) {
	  //   _CpTemp[i]= _a1cp1[i]*std::pow(temp,4)+_a2cp1[i]*std::pow(temp,3)+_a3cp1[i]*std::pow(temp,2)+_a4cp1[i]*temp+_a5cp1[i];
	  //}
	  //else {
	  //     _CpTemp[i]= _a1cp2[i]*std::pow(temp,4)+_a2cp2[i]*std::pow(temp,3)+_a3cp2[i]*std::pow(temp,2)+_a4cp2[i]*temp+_a5cp2[i];      
	  //}

	  
	  _CpTemp[i] /= _mmasses[i];
	  
	  
	  CpP += _CpTemp[i]*_ys[i];
      
     
      
      
   }
   
 
   CvP = CpP - (_Rgas * _massTtmTt);
   //It was
   // CvP = CpP - _Rgas;
   
   
   gamma = CpP / CvP ; 
   
   RoverMtot = _Rgas * _massTtmTt;
   
   //CFdouble gamma_old = 0.0;
   //gamma_old = 1 + (RoverMtot / CvP);
     
   soundSpeed = sqrt(gamma * temp * RoverMtot);
   
   //Old Sound Speed
   //CFdouble soundSpeed_old = 0.0; 
   //soundSpeed_old = sqrt(gamma*pressure/rho);
   
   //CFout << "frozengamma _mmasses O2 " << _mmasses[2] << "\n";
   //CFout << "frozengamma _Rgas " << _Rgas << "\n";
   //CFout << "frozengamma _massTtmTt " << _massTtmTt << "\n";
   //CFout << "frozengamma _Rgas * _massTtmTt " << _Rgas *_massTtmTt << "\n";
   //CFout << "frozengamma _CpTemp 0 " << _CpTemp[0] << "\n";
   //CFout << "frozengamma _CpTemp 1 " << _CpTemp[1] << "\n";
   //CFout << "frozengamma _CpTemp 2 " << _CpTemp[2] << "\n";
   //CFout << "frozengamma _CpTemp 3 " << _CpTemp[3] << "\n";
   //CFout << "frozengamma _CpTemp 4 " << _CpTemp[4] << "\n";
   //CFout << "frozengamma CpP " << CpP << "\n";
   //CFout << "frozengamma CvP " << CvP << "\n";
   //CFout << "frozengamma gamma " << gamma << "\n";
   //CFout << "frozengamma gamma_old " << gamma_old << "\n";
   //CFout << "frozengamma soundSpeed " << soundSpeed << "\n";
   //CFout << "frozengamma soundSpeed_old " << soundSpeed_old << "\n";
   
   
//    CFout << "gamma" << gamma << "\n";
//    CFout << "soundSpeed" << soundSpeed << "\n";
   
//CFout << "End of frozenGammaAndSoundSpeed" << "\n";
   
}

//////////////////////////////////////////////////////////////////////////////

void CombustionModelLibrary::setDensityEnthalpyEnergy(CFdouble& temp,
                                                //RealVector& tVec,
                                                CFdouble& pressure,
                                                RealVector& dhe)//,
                                                //bool storeExtraData)
{
  CFdouble hTot = 0.0;
  CFdouble _hsTot = 0.0;
  CFdouble RoverMtot = 0.0;
  CFint i;
  _massTtmTt = 0.0;
  _massTtmS = 0.0;
  
  //if NEW is used, uncomment
  //CFdouble _hformTot = 0.0;
  
   //CFout << "SetDensityEnthEnergy Start" << "\n";
//   CFout <<"Start of setDensityEnthalpyEnergy -- > dhe: " << "\n" << dhe << "\n";
//   CFdouble CvP = 0.0;
//   CFdouble CpP = 0.0;
  
//   CFout <<"inside setDensityEnthalpyEnergy" << "\n";

  //_meanmass = 0.0;
  
 // CFout << "temp at setDensityEnthalpyEnergy" << temp << "\n";
 
 
  dhe[0] = 0.0;
  dhe[1] = 0.0;
  dhe[2] = 0.0;
 

 switch (_thermoID){

      case 0: //Represents thermodynamic properties obtained through interpolation with Burcat polynomials
              // which give Htot = DHform_298 + integral(Cp dT);
      {
 
      for(i=0;i< _NS; ++i)
      {
 	 
	  if (temp <= 1000) {
		   _hsTemp[i]= _Rgas*temp*((_a1h1[i]) + (_a2h1[i] * temp)/2 + (_a3h1[i] * std::pow(temp,2))/3 + (_a4h1[i] * std::pow(temp,3))/4 + (_a5h1[i] * std::pow(temp,4))/5 + (_a6h1[i]/temp)); 
		
		   
		// old version for CEAS-Thermobuild coeffs, in data_old folder
		//CFout << "i-th enthalpy " << _hsTemp[2] << "\n";
		//_hsTemp[i]= _a1h1[i]*std::pow(temp,5)+_a2h1[i]*std::pow(temp,4)+_a3h1[i]*std::pow(temp,3)+_a4h1[i]*std::pow(temp,2)+_a5h1[i]*temp+_a6h1[i];
	  }
	  else {
		  _hsTemp[i]= _Rgas*temp*((_a1h2[i]) + (_a2h2[i] * temp)/2 + (_a3h2[i] * std::pow(temp,2))/3 + (_a4h2[i] * std::pow(temp,3))/4 + (_a5h2[i] * std::pow(temp,4))/5 + (_a6h2[i]/temp));
		  
		  
		  // old version for CEAS-Thermobuild coeffs, in data_old folder
	      //_hsTemp[i]= _a1h2[i]*std::pow(temp,5)+_a2h2[i]*std::pow(temp,4)+_a3h2[i]*std::pow(temp,3)+_a4h2[i]*std::pow(temp,2)+_a5h2[i]*temp+_a6h2[i];;
	  } 
//       CFout << "Before EnthPoly last line." << "\n";

      // Conversion from species' total enthalpy "per mole unit" to "per mass unit"
      //CFout << "hstemp " << _hsTemp << "\n";
	  _hsTemp[i] /= _mmasses[i];
	  
	  //OLD
	  _hsTot += _hsTemp[i] * _ys[i];
	  
	  
	  //NEW
	  //_hsTemp[i] = _hsTemp[i] - _hform[i]; //subtracted hform
	  
	  //_hsTot += _hsTemp[i]; //it was _hsTot += _hsTemp[i]*_ys[i];
	  
	  //_hformTot += _hform[i]*_ys[i]; //added
	  
	  //CFout << "setDensityEnthEnergy species enthalpy _htot " << hTot << "\n";
	  
	

      }
	
	  //CFout << "setdensity _ys " << _ys << "\n";
	  //CFout << "setdensity _hstemp " << _hsTemp << "\n";
	  //CFout << "setdensity _hform " << _hform << "\n";
	  
	  //CFout << "setdensity _hformTot " << _hformTot << "\n";
	  //CFout << "setdensity _hsTot " << _hsTot << "\n";
	  
	  //OLD
	  hTot = _hsTot;
	  
	  //NEW
	  //hTot = _hsTot + _hformTot; //added _hformTot
	  
      //CFout << "setdensity hTot " << hTot << "\n";
     
	  //CFout << "_Rgas " << _Rgas << "\n";
	  //CFout << "_mmasses " << _mmasses << "\n";
	  
	  
	  //CFout << "setDensityEnthEnergy species enthalpy _htot " << hTot << "\n";
      
      //CFout << "hTot" << hTot << "\n";
	  
	  for(i=0;i< _NS; ++i)
      {
		  _massTtmS[i] = _ys[i]/_mmasses[i];
		  _massTtmTt += _massTtmS[i];
		  
		  //CFout << "i = " << i << "\n";
		  //CFout << "_massTtmS 3 " << _massTtmS[3] << "\n";
		  //CFout << "_massTtmTt " << _massTtmTt << "\n";
	  }
      
      RoverMtot = _Rgas*_massTtmTt;
      
      dhe[0] = pressure/(RoverMtot*temp);
      
      dhe[1] = hTot;
      
      dhe[2] = dhe[1] - (pressure / dhe[0]);


	  //CFout << "RoverMtot " << RoverMtot << "\n";
	  
	  
//       CFout << "setDensityEnthalpyEnergy temperature" << temp <<"\n";
       //CFout <<"setDensityEnthalpyEnergy -- > density: " << dhe[0] << "\n";
//       CFout <<"setDensityEnthalpyEnergy -- > enthalpy: " << dhe[1] << "\n";
//       CFout <<"setDensityEnthalpyEnergy -- > energy: " << dhe[2] << "\n";
      //CFout << "SetDensityEnthEnergy End" << "\n";
      break;
      }
      
 }
 
}
      
//   for ( i = 0; i < _NS; ++i) {
//     
//         _massTtmS[i] = _ys[i]/_mmasses[i];
// 	_massTtmTt+=_massTtmS[i];
//           
//        	if (temp <= 1000) {
// 	     _CpTemp[i]= _a1cp1[i]*std::pow(temp,4)+_a2cp1[i]*std::pow(temp,3)+_a3cp1[i]*std::pow(temp,2)+_a4cp1[i]*temp+_a5cp1[i];
// 	  }
// 	  else {
// 	       _CpTemp[i]= _a1cp2[i]*std::pow(temp,4)+_a2cp2[i]*std::pow(temp,3)+_a3cp2[i]*std::pow(temp,2)+_a4cp2[i]*temp+_a5cp2[i];      
// 	  }
// 	  
// 	  _CpTemp[i] /= _mmasses[i];
// 	 
// 	  CpP += _CpTemp[i]*_ys[i];
// 	  
// 	  CvP = CpP - _Rgas;

	//CFout << "Mean Mass calculation Start" << "\n";
	
// 	_meanmass = 0.0;
// 	
// 	_meanmass += _mmasses[i];
// 	CFout << "Mean Mass:" << _meanmass << "\n";
// 	_meanmass = _meanmass/_NS;
	//CFout << "Mean Mass/_NS:" << _meanmass << "\n";
	
// 	_massTtmS[i]=_ys[i]/_mmasses[i];
//         _massTtmTt+=_massTtmS[i];
// 	
// 	_xs[i] = _ys[i]/_mmasses[i];
// 	_xs[i] = _xs[i]*(1/_massTtmTt);

//  CFdouble massTot=0.0;
// 
//   for (CFint is = 0; is < _NS; ++is) {
//     cf_always_assert (_ys[is] > 1.1);
// 
//     const CFreal mm = _ys[is]/_mmasses[is];
//     massTot += mm;
//     _xs[is] = mm;
//   }
//   _xs *= 1./massTot;
//       
// 	
//         //CFout << "_CpTemp - DONE" << "\n";
// 
// 	
//       CpP += _CpTemp[i]*_xs[i];
      


  //if (!storeExtraData) {
//     dhe[0] = (pressure*_massTtmTt)/(_Rgas*temp); //corrected from pressure/(_Rgas*_massTtmTt*temp);
//     dhe[1] = 0.0;
//     dhe[2] = 0.0;
    //dhe[3] = 0.0;
    
    
//     CFout <<"Start of setDensityEnthalpyEnergy -- > dhe: " << dhe[0] << "\n";
//     CFout <<"Start of setDensityEnthalpyEnergy -- > dhe: " << dhe[1] << "\n";
//     CFout <<"Start of setDensityEnthalpyEnergy -- > dhe: " << dhe[2] << "\n";
//     //CFout <<"Start of setDensityEnthalpyEnergy -- > dhe: " << dhe[3] << "\n";
//     CFout <<"Start of setDensityEnthalpyEnergy -- > dhe: " << dhe << "\n";

//     switch (_thermoID){
// 
//       case 0:
//       {
// 
//         for( i = 0; i < _NS; ++i) {
//           //dhe[3] += _ys[i]*energyVibSpeciesOverR(tVec[0],i);
//           hformTot += _ys[i]*_hform[i];
//         }
//         dhe[2] = ((CvP/_Rgas)+1)*temp; //+dhe[3];
//         dhe[2] *= _Rgas;
//         dhe[2] += hformTot;
//         //dhe[3] *= _Rgas;
//         dhe[1] = dhe[2] + pressure/dhe[0];
// 	
//    CFout << "End of setDensityEnthalpyEnergy." << "\n";
//     CFout <<"Start of setDensityEnthalpyEnergy -- > dhe0: " << dhe[0] << "\n";
//     CFout <<"Start of setDensityEnthalpyEnergy -- > dhe1: " << dhe[1] << "\n";
//     CFout <<"Start of setDensityEnthalpyEnergy -- > dhe2: " << dhe[2] << "\n";
//     //CFout <<"Start of setDensityEnthalpyEnergy -- > dhe3: " << dhe[3] << "\n";
//     CFout <<"hformTot: " << hformTot << "\n";
//     CFout <<"setDensityEnthalpyEnergy -- > dhe: " << dhe << "\n";
// 	
//       break;
//       }
  
//       }
// 
//   //}
// //   CFout <<"Start of setDensityEnthalpyEnergy -- > dhe: " << "\n" << dhe << "\n";
   //CFout << "SetDensityEnthEnergy End" << "\n";
//   
// }
  

      



//////////////////////////////////////////////////////////////////////////////

void CombustionModelLibrary::omegaContribution(CFdouble& temperature,
			             RealVector& tVec,
				     CFdouble& pressure,
                                     CFdouble& rho,
                                     const RealVector& ys,
				     const RealVector& mmasses,
                                     RealVector& omega)
{
  
 //omega = 0.0; //source term from chemistry
 CFint i;
//  CFout <<"omegaContribution start" << "\n";
 
 
 _Rgascal = 1.9858775; //Universal gas constant in cal/(mol K)

// //  for(i=0;i<_NS;i++){
// //    _ys[i]=ys[i];
// //  }

//  RealVector xs;

//   for (CFint is = 0; is < _NS; ++is) {
//     _ys[is] = ys[is];
// 
//     if (_ys[is] < 0.0) _ys[is] = 0.0;
//     cf_assert(_ys[is] <= 1.0)
//   }
//   
//     CFdouble massTot=0.0;
// 
//   for (CFint is = 0; is < _NS; ++is) {
//     cf_always_assert (ys[is] > 1.1);
// 
//     const CFreal mm = ys[is]/_mmasses[is];
//     massTot += mm;
//     xs[is] = mm;
//   }
//   xs *= 1./massTot;


 //NEED TO USE _ys or ys from ATDModel?????
 switch(_chemID)
 { 
   //CFout << "Temperature in omegaContribution" << temperature << "\n";
   case 0:
   {
     //Venkateswaran chemical model
     // Forward and backward reaction rates for Venkateswaran chemical model
      
      // k's are in the form k = A*exp(Ea/RT)
      // measurement units are;
      //A = [cm^3/(mol*s)]
      //Ea = [cal/mol]
      //_Rgascal = cal/(mol*K)
      // temperature = [K]
      //k = [cm^3/(mol*s)]
      
      _kfw1 = 5.9874*std::pow(10.,4)*exp(18001/(_Rgascal*temperature));
      
      _kfw2 = 8.4*std::pow(10.,11)*exp(63055/(_Rgascal*temperature)); //from Ibragimova et al. (ref [18])
      _kbw2 = 4.5*std::pow(10.,10)*exp(67116/(_Rgascal*temperature)); //from Ibragimova et al. (ref [18])
      
      // Conversion to k = [m^3/(kg*s)] for _kfw2 and _kbw2
      _kfw2 *= 1*std::pow(10.,-6);
      _kbw2 *= 1*std::pow(10.,-6);
      
      
      
     //Order of chemical species MATTERS (check respective .cmix file) 
     //C4H6
     _omegaTempDot[0] = -mmasses[0]*_kfw1*(rho*(ys[0]/mmasses[0]))*(rho*(ys[2]/mmasses[2]));
//      CFout << "omegadot [0] calculated" << "\n";
//      CFout << "omegaTempDot[0]"<< _omegaTempDot[0] << "\n";
     
     //H2O
     _omegaTempDot[1] = 3*mmasses[1]*_kfw1*(rho*(ys[0]/mmasses[0]))*(rho*(ys[2]/mmasses[2]));
//      CFout << "omegadot [1] calculated" << "\n";
//      CFout << "omegaTempDot[1]"<< _omegaTempDot[1] << "\n";
     
     //O2
     _omegaTempDot[2] = mmasses[2]*(-3.5*_kfw1*(rho*(ys[0]/mmasses[0]))*(rho*(ys[2]/mmasses[2]))-0.5*(_kfw2*(rho*(ys[3]/mmasses[3]))*sqrt((rho*(ys[2]/mmasses[2])))-_kbw2*(rho*(ys[4]/mmasses[4]))));
//      CFout << "omegadot [2] calculated" << "\n";
//      CFout << "omegaTempDot[2]"<< _omegaTempDot[2] << "\n";
     
     //CO
     _omegaTempDot[3] = mmasses[3]*(4*_kfw1*(rho*(ys[0]/mmasses[0]))*(rho*(ys[2]/mmasses[2]))-(_kfw2*(rho*(ys[3]/mmasses[3]))*sqrt((rho*(ys[2]/mmasses[2])))-_kbw2*(rho*(ys[4]/mmasses[4]))));
     
//      CFout << "omegadot [3] calculated" << "\n";
//      CFout << "omegaTempDot[3]"<< _omegaTempDot[3] << "\n";
     
     //CO2
     _omegaTempDot[4] = mmasses[4]*(_kfw2*(rho*(ys[3]/mmasses[3]))*sqrt((rho*(ys[2]/mmasses[2])))-_kbw2*(rho*(ys[4]/mmasses[4])));
     
//      CFout << "omegadot [4] calculated" << "\n";
//      CFout << "omegaTempDot[4]"<< _omegaTempDot[4] << "\n";
     
     
     for(i=0; i<_NS;i++)
     {
       omega[i] = _omegaTempDot[i];
     }
     break;
   }
   case 1:
   {
      //Jones-Lindstedt, 4 reaction 6 species model
      // Reaction rates for Jones4 model
      // all from Gariani PhD thesis
      _ka = 3.48*std::pow(10.,11)*exp(-30600/(_Rgascal*temperature));
      _kb = 9.11*std::pow(10.,13)*exp(-31200/(_Rgascal*temperature)); 
      _kc = 2.9*std::pow(10.,12)*exp(-19100/(_Rgascal*temperature));
      _kd = 2.8*std::pow(10.,18)*std::pow(temperature,-1)*exp(-43100/(_Rgascal*temperature));
//    Unused for Jones 4r 6s
//    _ke = 1.5*std::pow(10,9)*exp(-113000/(_Rgascal*temperatureerature));
//    _kf = 2.3*std::pow(10,22)*std::pow(temperatureerature,-3)*exp(-120000/(_Rgascal*temperatureerature));
      
      _Ra = _ka*((rho*(_ys[0]/mmasses[0]))*(rho*(_ys[1]/mmasses[1])));
      _Rb = _kb*(sqrt((rho*(_ys[0]/mmasses[0])))*std::pow((rho*(_ys[2]/mmasses[2])),1.25));
      _Rc = _kc*((rho*(_ys[4]/mmasses[4]))*(rho*(_ys[1]/mmasses[1])));
      _Rd = _kd*((std::pow(rho*(_ys[3]/mmasses[3]),0.25))*(std::pow(rho*(_ys[4]/mmasses[4]),1.5)));
//    Unused for Jones 4r 6s      
//    _Re = _ke*(mmasses[4]*_kfa*(rho*(_ys[4]/mmasses[4])));
//    2_Rf = _kf*(mmasses[4]*_kfa*(rho*(_ys[1]/mmasses[1])));*/
      
     //Order of chemical species MATTERS (check respective .cmix file) 
     //C4H6
     _omegaTempDot[0] = -mmasses[0]*(_Ra + _Rb);
     //H2O
     _omegaTempDot[1] =  mmasses[1]*(4*_Ra + _Rc - _Rd);
     //O2
     _omegaTempDot[2] =  mmasses[2]*(2*_Ra + 0.5*_Rd);
     //H2
     _omegaTempDot[3] =  mmasses[3]*(-3*_Ra - 7*_Rb - _Rc + _Rd);
     //CO
     _omegaTempDot[4] =  mmasses[4]*(-4*_Ra - 4*_Rb + _Rc);
     //CO2
     _omegaTempDot[5] =  mmasses[5]*(-_Rc);
     
      for(i=0; i<_NS;i++)
     {
       omega[i] = _omegaTempDot[i];
     }
     break;
   }
     
    case 2:
    {  
      //Jones-Lindstedt, 6 reaction 9 species model
      // Reaction rates for Jones4 model
      //all from Gariani PhD thesis
      _ka = 3.48*std::pow(10.,11)*exp(-30600/(_Rgascal*temperature));
      _kb = 9.11*std::pow(10.,13)*exp(-31200/(_Rgascal*temperature)); 
      _kc = 2.9*std::pow(10.,12)*exp(-19100/(_Rgascal*temperature));
      _kd = 2.8*std::pow(10.,18)*std::pow(temperature,-1)*exp(-43100/(_Rgascal*temperature));
      _ke = 1.5*std::pow(10.,9)*exp(-113000/(_Rgascal*temperature));
      _kf = 2.3*std::pow(10.,22)*std::pow(temperature,-3)*exp(-120000/(_Rgascal*temperature));
            
      _Ra = _ka*((rho*(ys[0]/mmasses[0]))*(rho*(ys[1]/mmasses[1])));
      _Rb = _kb*(sqrt((rho*(ys[0]/mmasses[0])))*std::pow((rho*(ys[2]/mmasses[2])),1.25));
      _Rc = _kc*((rho*(ys[4]/mmasses[4]))*(rho*(ys[1]/mmasses[1])));
      _Rd = _kd*((std::pow(rho*(ys[3]/mmasses[3]),0.25))*(std::pow(rho*(ys[4]/mmasses[4]),1.5)));
      _Re = _ke*(rho*(ys[4]/mmasses[4]));
      _Rf = _kf*(rho*(ys[1]/mmasses[1]));
            
     //Order of chemical species MATTERS (check respective .cmix file)
     //C4H6
     _omegaTempDot[0] = -mmasses[0]*(_Ra + _Rb);
     //H2O
     _omegaTempDot[1] =  mmasses[1]*(4*_Ra + _Rc - _Rd + _Rf);
     //H
     _omegaTempDot[2] = mmasses[2]*(_Rf);
     //O
     _omegaTempDot[3] = mmasses[3]*(2*_Re);
     //CO
     _omegaTempDot[4] =  mmasses[4]*(-4*_Ra - 4*_Rb + _Rc);
     //CO2
     _omegaTempDot[5] =  mmasses[5]*(-_Rc);
     //OH
     _omegaTempDot[6] =  mmasses[6]*(-_Rf);
     //O2
     _omegaTempDot[7] =  mmasses[7]*(2*_Ra + 0.5*_Rd + _Re);
     //H2
     _omegaTempDot[8] =  mmasses[8]*(-3*_Ra - 7*_Rb - _Rc + _Rd);
     
      for(i=0; i<_NS;i++)
     {
       omega[i] = _omegaTempDot[i];
     }
     break;
    }
     
 }
 
//  CFout << "End of omegaContribution." << "\n";
 
//return omegadot;

}

//////////////////////////////////////////////////////////////////////////////

void CombustionModelLibrary::getMassProductionTerm(CFdouble& temp, 
						   RealVector& tVec, 
						   CFdouble& pressure, 
						   CFdouble& rho, 
						   const RealVector& ys,
						   bool flagJac, 
						   RealVector& omega, 
						   RealMatrix& jacobian)
{
  //CFdouble omegadot;
  
// OLD VERSION JUST CALLING omegaContribution function
//   CFout << "getMassProductionTerm Begins" << "\n";
//   omega = 0.0;
//   
//   omegaContribution(temperature,rho,ys,_mmasses,omega);


// NEW VERSION - FULL omegaContribution re-written

 CFint i;
//  CFout <<"getMassProductionTerm starts" << "\n";
 
 _Rgascal = 1.9858775; //Universal gas constant in cal/(mol*K)

 //NEED TO USE _ys or ys from ATDModel?????
 //CFout << "Temperature in getMassProductionTerm" << temp << "\n";
 switch(_chemID)
 { 
   case 0:
    // CFout << "Temperature in getMassProductionTerm" << temp << "\n";
   {//Venkateswaran chemical model
      // Forward and backward reaction rates for Venkateswaran chemical model
      
//       CFout << "inside getomegadot function - case 0" << "\n";
      
      // k's are in the form k = A*exp(Ea/RT)
      // measurement units are;
      //A = [cm^3/(mol*s)]
      //Ea = [cal/mol]
      //_Rgascal = cal/(mol*K)
      // temperature = [K]
      //k = [cm^3/(mol*s)]
	  //_mmasses[i] = [kg/mol] (here, converted after reading)
	  //rho[i] = [kg/m^3]
      
 ////////////////////////////////////////////////////// OLD CHEMISTRY//////////////////////////////////////////
	  // OLD value
      //_kfw1 = 5.9874*std::pow(10,4)*exp(18001/(_Rgascal*temp));
	  //CFdouble _Ad = 2.904884967 * std::pow(10.,13);  // *N for ignition when restart from developed fluid field
	  //CFdouble _Ed = 51000;	//cal/mol
	  //CFdouble _As = 5.987414172 * std::pow(10.,4);   // *N for ignition when restart from developed fluid field
	  //CFdouble _Es = 18000;	//cal/mol
	  
	  //CFdouble _kd = _Ad * exp(-_Ed/(_Rgascal*temp));
	  //CFdouble _ks = _As * exp(-_Es/(_Rgascal*temp));
	  
	  
	   //_kfw1 = ((_ks*(pressure/(_Rgas*298))) + (_kd * (rho*ys[0]/_mmasses[0])))/((rho*ys[0]/_mmasses[0]) * (1 + (_kd/_ks)) + (pressure/(_Rgas*298))); 
//	   //multiply _kfw1 *N for ignton when restart from developed fluid field
//	   //_kfw1 = _kfw1*10;

//           //tryout: uncomment to kill first reaction (cfr with OpenFoam)
//           //_kfw1 = 0;
//          
//	   
//	   
//	  // OLD values
//      //_kfw2 = 8.4*std::pow(10,11)*exp(63055/(_Rgascal*temp)); //from Ibragimova et al. (ref [18])
//      //_kbw2 = 4.5*std::pow(10,10)*exp(67116/(_Rgascal*temp)); //from Ibragimova et al. (ref [18])
//	   
//	  
//	  // NEW values from Westbrook and Dryer, 1981
//	  CFdouble _concC4H6 = (rho*(ys[0]/_mmasses[0]));
//	  CFdouble _concH2O = (rho*(ys[1]/_mmasses[1]));
//	  CFdouble _concO2 = (rho*(ys[2]/_mmasses[2]));
//	  CFdouble _concCO = (rho*(ys[3]/_mmasses[3]));
//	  CFdouble _concCO2 = (rho*(ys[4]/_mmasses[4]));
//	  
//	  CFdouble _A2 = 1 * std::pow(10,14.6); // *N for ignition when restart from developed fluid field
//	  CFdouble _A3 = 5 * std::pow(10,8);  // *N for ignition when restart from developed fluid field
//	  //old version, not correct, from Westbrook-Dryer, 1981
//	  //pre exponential A (8.8e11 is for C4H10)
//	  //CFdouble _A3 = 8.8 * std::pow(10,11);  // *N for ignition when restart from developed fluid field
//	  
//	  CFdouble _Ea2 = -40000;
//				  
//	  _kfw2 = _A2 * exp(_Ea2/(_Rgascal*temp)) * _concCO * std::pow(_concH2O,0.5) * std::pow(_concO2,0.25);
//	  _kbw2 = _A3 * exp(_Ea2/(_Rgascal*temp)) * _concCO2; 
//	  
//	   
//      // Conversion to k = [m^3/(kg*s)] for _kfw2 and _kbw2
//      _kfw2 *= 1*std::pow(10,-6);
//      _kbw2 *= 1*std::pow(10,-6);

//      //tryout: uncomment to kill second reaction
//      _kfw2 = 0;
//      _kbw2 = 0;
//      
//      
//      
////       CFout << "Reaction constants set" << "\n";
//      
////       CFout << "_kfw1" << _kfw1 <<"\n";
////       CFout << "_kfw2" << _kfw1 <<"\n";
////       CFout << "_kbw2" << _kfw1 <<"\n";
////       
////       //ATDModelLibrary::setSpeciesFractions(ys);
////       
////       CFout << "ys[0]" << ys[0] <<"\n";
////       CFout << "ys[1]" << ys[1] <<"\n";
////       CFout << "ys[2]" << ys[2] <<"\n";
////       CFout << "ys[3]" << ys[3] <<"\n";
////       CFout << "ys[4]" << ys[4] <<"\n";
//      
//      
//      
//      
//     //Order of chemical species MATTERS (check respective .cmix file) 
//     //C4H6
//     //_omegaTempDot[0] = -_mmasses[0]*_kfw1*(rho*(ys[0]/_mmasses[0]))*(rho*(ys[2]/_mmasses[2]));
//	 _omegaTempDot[0] = -_mmasses[0]*_kfw1*_concC4H6 * _concO2;
////      CFout << "omegadot [0] calculated" << "\n";
////      CFout << "omegaTempDot[0]"<< _omegaTempDot[0] << "\n";
//     
//     //H2O
//     //_omegaTempDot[1] = 3*_mmasses[1]*_kfw1*(rho*(ys[0]/_mmasses[0]))*(rho*(ys[2]/_mmasses[2]));
//	 _omegaTempDot[1] = 3*_mmasses[1]*_kfw1*_concC4H6 *_concO2;
////      CFout << "omegadot [1] calculated" << "\n";
////      CFout << "omegaTempDot[1]"<< _omegaTempDot[1] << "\n";
//     
//     //O2
//     //_omegaTempDot[2] = _mmasses[2]*(-3.5*_kfw1*(rho*(ys[0]/_mmasses[0]))*(rho*(ys[2]/_mmasses[2]))-0.5*(_kfw2*(rho*(ys[3]/_mmasses[3]))*sqrt((rho*(ys[2]/_mmasses[2])))-_kbw2*(rho*(ys[4]/_mmasses[4]))));
//	 _omegaTempDot[2] = _mmasses[2]*((-3.5*_kfw1*_concC4H6 * _concO2) - 0.5*(_kfw2*_concCO * sqrt(_concO2)-_kbw2*_concCO2));
//	 
////      CFout << "omegadot [2] calculated" << "\n";
////      CFout << "omegaTempDot[2]"<< _omegaTempDot[2] << "\n";
//     
//     //CO
//     //_omegaTempDot[3] = _mmasses[3]*(4*_kfw1*(rho*(ys[0]/_mmasses[0]))*(rho*(ys[2]/_mmasses[2]))-(_kfw2*(rho*(ys[3]/_mmasses[3]))*sqrt((rho*(ys[2]/_mmasses[2])))-_kbw2*(rho*(ys[4]/_mmasses[4]))));
//	   _omegaTempDot[3] = _mmasses[3] * ((4 * _kfw1 * _concC4H6 * _concO2) - (_kfw2 * _concCO * sqrt(_concO2) - _kbw2 * _concCO2));
//     
////      CFout << "omegadot [3] calculated" << "\n";
////      CFout << "omegaTempDot[3]"<< _omegaTempDot[3] << "\n";
//     
//     //CO2
//     //_omegaTempDot[4] = _mmasses[4]*(_kfw2*(rho*(ys[3]/_mmasses[3]))*sqrt((rho*(ys[2]/_mmasses[2])))-_kbw2*(rho*(ys[4]/_mmasses[4])));
//	   _omegaTempDot[4] = _mmasses[4]*((_kfw2 * _concCO * sqrt(_concO2)) - _kbw2 * _concCO2);
//     
////      CFout << "omegadot [4] calculated" << "\n";
////      CFout << "omegaTempDot[4]"<< _omegaTempDot[4] << "\n";
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////     
     
	 //Pre-exponential and activation energy for C4H6 oxidation reaction, from Westbrook-Dryer 1981
	 // only forward
	 CFdouble _Af1 = 8.8 * std::pow(10.,11);
	 CFdouble _Eaf1 = 30000;
         
	 //MODDING ARTIFICIALLY _Af1 - COMMENT/UNCOMMENT AS NEEDED
         _Af1 = 15 * _Af1;
	 
	 //Pre-exponential and activation energy for CO oxidation reaction, from Westbrook-Dryer 1981
	 // forward
	 CFdouble _Af2 = 1 * std::pow(10.,14.6);
	 CFdouble _Eaf2 = 40000;
	 // backward  
	 CFdouble _Ab2 = 5 * std::pow(10.,8);
	 CFdouble _Eab2 = 40000;
	 
	 //Concentrations in mol/m^3
	 CFdouble _concC4H6 = (rho*(ys[0]/_mmasses[0]));
	 CFdouble _concH2O = (rho*(ys[1]/_mmasses[1]));
	 CFdouble _concO2 = (rho*(ys[2]/_mmasses[2]));
	 CFdouble _concCO = (rho*(ys[3]/_mmasses[3]));
	 CFdouble _concCO2 = (rho*(ys[4]/_mmasses[4]));
	  
	 // Concentrations in mol/cm^3 (NEEDED for use in all reaction _k)
	 CFdouble _concC4H6cm = _concC4H6 * std::pow(10.,-6);
	 CFdouble _concH2Ocm = _concH2O * std::pow(10.,-6);
	 CFdouble _concO2cm = _concO2 * std::pow(10.,-6);
	 CFdouble _concCOcm = _concCO * std::pow(10.,-6);
	 CFdouble _concCO2cm = _concCO2 * std::pow(10.,-6);
	 
	 // Reaction rates in mol/(cm^3 s)
	 _kfw1 = _Af1 * exp(-_Eaf1/(_Rgascal*temp)) * std::pow(_concC4H6cm,0.15) * std::pow(_concO2cm,1.6);
	 
	 _kfw2 = _Af2 * exp(-_Eaf2/(_Rgascal*temp)) * _concCOcm * std::pow(_concH2Ocm,0.5) * std::pow(_concO2cm,0.25);
	 _kbw2 = _Ab2 * exp(-_Eab2/(_Rgascal*temp)) * _concCO2cm; 
	 
	 // Converting reaction rates back to SI units (mol/(m^3 s))
	 _kfw1 = _kfw1 * std::pow(10.,6);
	 
	 _kfw2 = _kfw2 * std::pow(10.,6); 
	 _kbw2 = _kbw2 * std::pow(10.,6); 
	 
	 
	 //Order of chemical species MATTERS (check respective .cmix file) 
     //C4H6
	 _omegaTempDot[0] = -_mmasses[0] * _kfw1;
     
     //H2O
	 _omegaTempDot[1] = 3 * _mmasses[1] * _kfw1;
	 
	 //O2
	 _omegaTempDot[2] = _mmasses[2] * (-3.5 * _kfw1 - 0.5 * _kfw2 + 0.5 * _kbw2);
	 
	 //CO
	 _omegaTempDot[3] = _mmasses[3] * (4 * _kfw1 - _kfw2 + _kbw2);
	 
	 //CO2
	 _omegaTempDot[4] = _mmasses[4] * (_kfw2 - _kbw2);
		 
	   
     for(i=0; i<_NS;i++)
     {
       omega[i] = _omegaTempDot[i];
        //CFout << "omegadot vector" << omega[i] << "\n";
     }
	 
	// ANALYTICAL JACOBIAN FOR OLD CHEMISTRY - NO MORE CORRECT
	// flagJac == false;
	// if (flagJac==true)
	// {
	//	 const CFuint dim = PhysicalModelStack::getActive()->getDim();
	//	 const CFuint nbEqs = jacobian.nbRows();
	//	 assert(dim == nbEqs - _NS - 1);
	//	 const CFuint sizeJ = nbEqs - dim;
	//	 const CFuint startT = _NS + dim;
	//	 
	//	 jacobian = 0.0;
	//	 
	//// Assembling Jacobian matrix of the source term corresponding to
    //// primitive variables: rho_i, u, v,( w,) T, (Tv, Te), (k, omega)

    //// suppose N = _NS-1 =>
    //// jacobian (nbEq*nbEq) must contain the following derivatives in this order:

    //// dOm(0)/dRho_s dOm(0)/dV dOm(0)/dT dOm(0)/dTv_m dOm(0)/dTe dOm(0)/dk dOm(0)/domega
    //// ...
    //// dOm(N)/dRho_s dOm(N)/dV dOm(N)/dT dOm(N)/dTv_m dOm(N)/dTe dOm(N)/dk dOm(0)/domega

    //// -----_NS----- ---DIM---   ----1----  --1 for k---   -- 1 for omega---
    //// with        dOm(s)/dV = 0            Om(s)/dk = 0   dOm(s)/domega = 0
	//	 
	//	 CFdouble _dkfw1_drho0 = (exp(_Ed / (_Rgas*temp))) / (  (_mmasses[0] * pressure)/(_Rgas * 298) + ((rho*ys[0]) * (((_Ad * exp(_Ed / (_Rgas * temp)) )/(_As * exp(_Es / (_Rgas * temp)))) + 1))  );
	//	 _dkfw1_drho0 += - ((((rho*ys[0])/_mmasses[0]) * exp(_Ed/(_Rgas * temp)) + (_As * pressure * exp(_Es/(_Rgas*temp)))/(_Rgas * 298)) *  (((_Ad * exp(_Ed / (_Rgas * temp)) )/(_As * exp(_Es / (_Rgas * temp)))) + 1)  )/ ((  (_mmasses[0] * pressure)/(_Rgas * 298) + ((rho*ys[0]) * (((_Ad * exp(_Ed / (_Rgas * temp)) )/(_As * exp(_Es / (_Rgas * temp)))) + 1))  ) * (  (_mmasses[0] * pressure)/(_Rgas * 298) + ((rho*ys[0]) * (((_Ad * exp(_Ed / (_Rgas * temp)) )/(_As * exp(_Es / (_Rgas * temp)))) + 1))  ));
	//					 
	//	 CFdouble _dkfw1_dT = - (_Ed * (rho*ys[0]) * exp(_Ed/(_Rgas*temp))) / (_mmasses[0] * _Rgas * temp); 
	//     _dkfw1_dT += - (_As * _Es * pressure * exp(_Es/(_Rgas*temp)) / (_Rgas * _Rgas * temp * temp * 298));
	//	 
	//	 CFdouble _dkfw2_dT = -_Ea2 / (_Rgas * temp * temp) * _kfw2;
	//	 CFdouble _dkbw2_dT = -_Ea2 / (_Rgas * temp * temp) * _kbw2;
	// 
	//	 //Non-zero derivatives for omega[0]
	//	 //jacobian(0,0) = - ((rho*ys[2]) / _mmasses[2]) * (_kfw1 + (rho*ys[0]) * _dkfw1_drho0);
	//	 //jacobian(2,0) = -_kfw1 * rho * ys[0] / _mmasses[2];				 
	//	 //jacobian(0,startT) = - ((rho * ys[0]) * (rho * ys[2]) / _mmasses[2]) * _dkfw1_dT;

	//	 //Correction of ys[1] and ys[2] to keep them > 0
	//	 //ys[1]= ys[1] + 1*std::pow(10,-16);
	//	 //ys[2]= ys[2] + 1*std::pow(10,-16);
	//	 
	//	 
    //     jacobian(0,0) = - _concO2 * (_kfw1 + (rho*ys[0]) * _dkfw1_drho0);
	//	 jacobian(2,0) = -_kfw1 * rho * ys[0] / _mmasses[2];
    //     jacobian(0,startT) = - ((rho * ys[0]) * _concO2) * _dkfw1_dT;

	//	 
	//	 //Non-zero derivatives for omega[1]
	//	 //jacobian(0,1) = 3 * ((rho*ys[2]) / _mmasses[2]) * (_kfw1 + (rho*ys[0]) * _dkfw1_drho0);
	//	 //jacobian(2,1) = 3 * _kfw1 * rho * ys[0] / _mmasses[2];				 
	//	 //jacobian(1,startT) = 3 * ((rho * ys[0]) * (rho * ys[2]) / _mmasses[2]) * _dkfw1_dT;
	//	 
	//	 jacobian(0,1) = 3 * (_concO2 * (_kfw1 + (rho*ys[0]) * _dkfw1_drho0));
	//	 jacobian(2,1) = 3 * _kfw1 * rho * ys[0] / _mmasses[2];				 
	//	 jacobian(1,startT) = 3 * ((rho * ys[0]) * _concO2) * _dkfw1_dT;
	//	 
	//	 
	//	 //Non-zero derivatives for omega[2]
	//	 jacobian(0,2) = - 3.5 * rho * (ys[2] / _mmasses[0]) * (rho * ys[0] * _dkfw1_drho0 + _kfw1);
	//     jacobian(1,2) = - (_mmasses[2] / (rho*ys[1])) * _kfw2 * _concCO * sqrt(_concO2);
    //     jacobian(2,2) = - 3.5 * _concC4H6 * _kfw1 - (3/8) * _concCO * _kfw2 * sqrt(1/((rho*(ys[2]/_mmasses[2]))));
    //     jacobian(3,2) = - (_mmasses[2] / _mmasses[3]) * sqrt(_concO2) * _kfw2;
    //     jacobian(4,2) = (_mmasses[2] / _mmasses[4]) * _kbw2;
    //     jacobian(2,startT) = _mmasses[2] * ( - 3.5 * (_concC4H6 * _concO2 * _dkfw1_dT) - 0.5 * ((_concCO * sqrt(_concO2) * _dkfw2_dT) - (_concCO2 * _dkbw2_dT)));


    //     //Non-zero derivatives for omega[3]
	//	 jacobian(0,3) = 4 * _mmasses[3] * _concO2 * (_concC4H6 * _dkfw1_drho0 + _kfw1/_mmasses[0]);
	//     jacobian(1,3) = -0.5 * _mmasses[3] * _concCO * sqrt(_concO2) * _kfw2/(rho*ys[1]);
    //     jacobian(2,3) = 4 * _concC4H6 * (_mmasses[3]/_mmasses[2]) * _kfw1 - (3/4) * (rho*ys[3]) * (1/(sqrt(_mmasses[2]))) * (1/(sqrt(rho*ys[2]))) * _kfw2;
    //     jacobian(3,3) = - 2 * sqrt(_concO2) * _kfw2;
    //     jacobian(4,3) = 2 * (_mmasses[3]/_mmasses[4]) * _kbw2;
    //     jacobian(3,startT) = _mmasses[3] * ( 4 * (_concC4H6 * _concO2 * _dkfw1_dT) - ((_concCO * sqrt(_concO2) * _dkfw2_dT) - (_concCO2 * _dkbw2_dT)));

    //         
    //     //Non-zero derivatives for omega[4]
    //     jacobian(1,4) = 0.5 * _mmasses[4] * _concCO * sqrt(_concO2) * _kfw2/(rho*ys[1]);
    //     jacobian(2,4) = (3/4) * _concCO * _kfw2 * _mmasses[4] * (1/(sqrt(_mmasses[2]))) * (1/(sqrt(rho*ys[2])));
    //     jacobian(3,4) = 2 * (_mmasses[4]/_mmasses[3]) * sqrt(_concO2) * _kfw2;
    //     jacobian(4,4) = -2 * _kbw2;
    //     jacobian(4,startT) = _mmasses[4] * ((_concCO * sqrt(_concO2) * _dkfw2_dT) - (_concCO2 * _dkbw2_dT));
		 
	 //}
          
     break;
   }
   case 1:
     //CFout << "Temperature in getMassProductionTerm" << temp << "\n";
   {
      //Jones-Lindstedt, 4 reaction 6 species model
      // Reaction rates for Jones4 model
      // all from Gariani PhD thesis
	  //Species order for Jones4-6
	  //  0  -  1  -  2 - 3 - 4 - 5
	  // C4H6  H2O   O2  H2  CO  CO2
      _ka = 3.48*std::pow(10.,11)*exp(-30600/(_Rgascal*temp));
      _kb = 9.11*std::pow(10.,13)*exp(-31200/(_Rgascal*temp)); 
      _kc = 2.9*std::pow(10.,12)*exp(-19100/(_Rgascal*temp));
      _kd = 2.8*std::pow(10.,18)*std::pow(temp,-1)*exp(-43100/(_Rgascal*temp));
//    Unused for Jones 4r 6s
//    _ke = 1.5*std::pow(10,9)*exp(-113000/(_Rgascal*temperature));
//    _kf = 2.3*std::pow(10,22)*std::pow(temperature,-3)*exp(-120000/(_Rgascal*temperature));
      
      _Ra = _ka*((rho*(_ys[0]/_mmasses[0]))*(rho*(_ys[1]/_mmasses[1])));
      _Rb = _kb*(sqrt((rho*(_ys[0]/_mmasses[0])))*std::pow((rho*(_ys[2]/_mmasses[2])),1.25));
      _Rc = _kc*((rho*(_ys[4]/_mmasses[4]))*(rho*(_ys[1]/_mmasses[1])));
      _Rd = _kd*((std::pow(rho*(_ys[3]/_mmasses[3]),0.25))*(std::pow(rho*(_ys[2]/_mmasses[2]),1.5)));
//    Unused for Jones 4r 6s      
//    _Re = _ke*(_mmasses[4]*_kfa*(rho*(_ys[4]/_mmasses[4])));
//    2_Rf = _kf*(_mmasses[4]*_kfa*(rho*(_ys[1]/_mmasses[1])));*/
      
     //Order of chemical species MATTERS (check respective .cmix file) 
     //C4H6
     _omegaTempDot[0] =  -_mmasses[0]*(_Ra + _Rb);
     //H2O
     _omegaTempDot[1] =  -_mmasses[1]*(4*_Ra + _Rc - _Rd);
     //O2
     _omegaTempDot[2] =  -_mmasses[2]*(2*_Rb + 0.5*_Rd);
     //H2
     _omegaTempDot[3] =  -_mmasses[3]*(- 7*_Ra -3*_Rb - _Rc + _Rd);
     //CO
     _omegaTempDot[4] =  -_mmasses[4]*(-4*_Ra - 4*_Rb + _Rc);
     //CO2
     _omegaTempDot[5] =  -_mmasses[5]*(-_Rc);
     
      for(i=0; i<_NS;i++)
     {
       omega[i] = _omegaTempDot[i];
     }
     break;
   }
     
    case 2:
      //CFout << "Temperature in getMassProductionTerm" << temp << "\n";
    {  
      //Jones-Lindstedt, 6 reaction 9 species model
      // Reaction rates for Jones4 model
      //all from Gariani PhD thesis
		//Species order for Jones4-6
		//  0  -  1  - 2 - 3 - 4 - 5  - 6 - 7 - 8
		// C4H6  H2O   H   O  CO  CO2  OH  O2  H2
      _ka = 3.48*std::pow(10.,11)*exp(-30600/(_Rgascal*temp));
      _kb = 9.11*std::pow(10.,13)*exp(-31200/(_Rgascal*temp)); 
      _kc = 2.9*std::pow(10.,12)*exp(-19100/(_Rgascal*temp));
      _kd = 2.8*std::pow(10.,18)*std::pow(temp,-1)*exp(-43100/(_Rgascal*temp));
      _ke = 1.5*std::pow(10.,9)*exp(-113000/(_Rgascal*temp));
      _kf = 2.3*std::pow(10.,22)*std::pow(temp,-3)*exp(-120000/(_Rgascal*temp));
            
      _Ra = _ka*((rho*(ys[0]/_mmasses[0]))*(rho*(ys[1]/_mmasses[1])));
      _Rb = _kb*(sqrt((rho*(ys[0]/_mmasses[0])))*std::pow((rho*(ys[2]/_mmasses[2])),1.25));
      _Rc = _kc*((rho*(ys[4]/_mmasses[4]))*(rho*(ys[1]/_mmasses[1])));
      _Rd = _kd*((std::pow(rho*(ys[8]/_mmasses[8]),0.25))*(std::pow(rho*(ys[7]/_mmasses[7]),1.5)));
      _Re = _ke*(rho*(ys[7]/_mmasses[7]));
      _Rf = _kf*(rho*(ys[1]/_mmasses[1]));
            
     //Order of chemical species MATTERS (check respective .cmix file)
     //C4H6
     _omegaTempDot[0] =  -_mmasses[0]*(_Ra + _Rb);
     //H2O
     _omegaTempDot[1] =  -_mmasses[1]*(4*_Ra + _Rc - _Rd + _Rf);
     //H
     _omegaTempDot[2] = -_mmasses[2]*(-_Rf);
     //O
     _omegaTempDot[3] = -_mmasses[3]*(2*-_Re);
     //CO
     _omegaTempDot[4] =  -_mmasses[4]*(-4*_Ra - 4*_Rb + _Rc);
     //CO2
     _omegaTempDot[5] =  -_mmasses[5]*(-_Rc);
     //OH
     _omegaTempDot[6] =  -_mmasses[6]*(-_Rf);
     //O2
     _omegaTempDot[7] =  -_mmasses[7]*(2*_Rb + 0.5*_Rd + _Re);
     //H2
     _omegaTempDot[8] =  -_mmasses[8]*(- 7*_Ra -3*_Rb - _Rc + _Rd);
     
      for(i=0; i<_NS;i++)
     {
       omega[i] = _omegaTempDot[i];
     }
     break;
    }
     
 }
 
// CFout <<"getMassProductionTerm ends:" << " omega" << omega << "\n";

  

// CFint i;
//   omega=0.0;
  
//   for(i=0;i<(CFint) _ChemReactArray.size();i++)
//     {
//       _ChemReactArray[i].omegaContribution(temperature,tVec,pressure,rho,ys,__mmasses,omega);
//     }
//   omega*=__mmasses;
  
//   if (flagJac)
//     {
//       CFout << "Function not implemented: CombustionModelLibrary::getMassProductionTerm(...,flagJac=true ,...) \n";
//       throw Common::NotImplementedException(FromHere(),"CombustionModelLibrary::getMassProductionTerm()");
//     }

}

//////////////////////////////////////////////////////////////////////////////




void CombustionModelLibrary::getRhoUdiff(CFdouble& temp,
					 CFdouble& pressure,
					 RealVector& normConcGradients,
					 CFreal* tVec,
					 RealVector& rhoUdiff,
					 bool fast)
{
//   CFout << "getRhoUdiff Begins" << "\n";
  CFdouble lambda = 0.0;
  CFdouble CpP = 0.0;
  CFint i = 0;
//   _massTtmTt = 0.0;
//   _massTtmS = 0.0;
//   //_meanmass = 0.0;
  

  
  
  
  //_D = 0.0;

  lambda = lambdaNEQ(temp,pressure);
  

  for ( i = 0; i < _NS; ++i) {
    
//      _massTtmS[i] = _ys[i]/_mmasses[i];
// 	_massTtmTt+=_massTtmS[i];

//SPECIFIC HEAT Cp CALCULATION
	  
	  if (temp <= 1000) {
		  _CpTemp[i] = _Rgas * (_a1cp1[i] + (_a2cp1[i] * temp) + (_a3cp1[i] * std::pow(temp,2)) + (_a4cp1[i] * std::pow(temp,3)) + (_a5cp1[i] * std::pow(temp,4))); 
		  
		  
	  }
	  else {
		  _CpTemp[i] = _Rgas * (_a1cp2[i] + (_a2cp2[i] * temp) + (_a3cp2[i] * std::pow(temp,2)) + (_a4cp2[i] * std::pow(temp,3)) + (_a5cp2[i] * std::pow(temp,4))); 
		  
	  }
	  
	  // old version for CEAS-Thermobuild coeffs, in data_old folder	
      // 	if (temp <= 1000) {
	  //   _CpTemp[i]= _a1cp1[i]*std::pow(temp,4)+_a2cp1[i]*std::pow(temp,3)+_a3cp1[i]*std::pow(temp,2)+_a4cp1[i]*temp+_a5cp1[i];
	  //}
	  //else {
	  //     _CpTemp[i]= _a1cp2[i]*std::pow(temp,4)+_a2cp2[i]*std::pow(temp,3)+_a3cp2[i]*std::pow(temp,2)+_a4cp2[i]*temp+_a5cp2[i];      
	  //}
      
	  // Conversion from "per mole unit" to "per mass unit" of Cp
	  _CpTemp[i] /= _mmasses[i];
	  
	  CpP += _CpTemp[i]*_ys[i];
  }
	  
  for ( i = 0; i < _NS; ++i) {
	  rhoUdiff[i]=-normConcGradients[i]*(lambda)/(CpP*_Le);
	  
	  //CFout << "getrhoudiff rhoudiff[i] " << rhoUdiff[i] << "\n";
	  //It was:
	  //rhoUdiff[i]=-normConcGradients[i]*(lambda*_Rgas)/(CpP*_Le); 
  }  
 
	  

//  CFdouble massTot=0.0;
// 
//   for (CFint is = 0; is < _NS; ++is) {
//     cf_always_assert (_ys[is] > 1.1);
// 
//     const CFreal mm = _ys[is]/_mmasses[is];
//     massTot += mm;
//     _xs[is] = mm;
//   }
//   _xs *= 1./massTot;
      
//   for ( i = 0; i < _NS; ++i) {    
//   CpP += _CpTemp[i]*_xs[i]; // correct to _ys
//   }

  //_D = (lambda/(rho*CpP*_Le));
  
    
//   for(i=0;i< _NS;i++)
//   {
//       //rhoUdiff[i]=-normConcGradients[i]*(lambda/(_Le*cpP));
//       rhoUdiff[i]=-normConcGradients[i]*(lambda*_Rgas)/(CpP*_Le); //remove _Rgas
//   }
// CFout << "normConcGradients: "<< normCon''cGradients << "\n";
// CFout << "lambda: "<< lambda << "\n";
// CFout << "CpP: "<< CpP << "\n";
// CFout << "getRhoUdiff Ends: "<< rhoUdiff << "\n";
}

//////////////////////////////////////////////////////////////////////////////

void CombustionModelLibrary::getDij_fick(RealVector& dx,
                                  CFdouble& pressure,
                                  CFdouble& temperature,
                                  CFreal& Diff_coeff,
                                  RealMatrix& Dij,
                                  RealVector& rhoUdiff)
{
  //not implemented in mutation 2.0
  // look at plugins/MutationI/MutationLibrary.cxx
  CFout << "Function not implemented: CombustionModelLibrary::getDij_fick() \n";
  throw Common::NotImplementedException(FromHere(),"CombustionModelLibrary::getDij_fick()");
}

//////////////////////////////////////////////////////////////////////////////

void CombustionModelLibrary::getSpeciesTotEnthalpies(CFdouble& temp,
                                               RealVector& tVec,
                                               CFdouble& pressure,
                                               RealVector& hsTot,
                                               RealVector* hsVib,
                                               RealVector* hsEl)



{
  //   CFout << "getSpeciesTotEnthalpies Begins" << "\n";
  switch(_thermoID)
  {
  case 0: //Polynomial interpolation of Cp - Represents thermodynamic properties obtained through interpolation with Burcat polynomials
          // which give Htot = DHform_298 + integral(Cp dT);

  CFint i;//,j;
  for(i=0;i< _NS;i++)
    {
// 	_massTtmS[i]=_ys[i]/_mmasses[i];
//         _massTtmTt+=_massTtmS[i];
// 	
		if (temp <= 1000) {
			_hsTemp[i]= _Rgas*temp*((_a1h1[i]) + (_a2h1[i] * temp)/2 + (_a3h1[i] * std::pow(temp,2))/3 + (_a4h1[i] * std::pow(temp,3))/4 + (_a5h1[i] * std::pow(temp,4))/5 + (_a6h1[i]/temp)); 
									// old version for CEAS-Thermobuild coeffs, in data_old folder
									//_hsTemp[i]= _a1h1[i]*std::pow(temp,5)+_a2h1[i]*std::pow(temp,4)+_a3h1[i]*std::pow(temp,3)+_a4h1[i]*std::pow(temp,2)+_a5h1[i]*temp+_a6h1[i];
		}
		else {
			_hsTemp[i]= _Rgas*temp*((_a1h2[i]) + (_a2h2[i] * temp)/2 + (_a3h2[i] * std::pow(temp,2))/3 + (_a4h2[i] * std::pow(temp,3))/4 + (_a5h2[i] * std::pow(temp,4))/5 + (_a6h2[i]/temp));
									// old version for CEAS-Thermobuild coeffs, in data_old folder
									//_hsTemp[i]= _a1h2[i]*std::pow(temp,5)+_a2h2[i]*std::pow(temp,4)+_a3h2[i]*std::pow(temp,3)+_a4h2[i]*std::pow(temp,2)+_a5h2[i]*temp+_a6h2[i];;
		} 
		
//       CFout << "Before EnthPoly last line." << "\n";

      // Conversion from species' total enthalpy "per mole unit" to "per mass unit"
      _hsTemp[i] /= _mmasses[i];
	  
      
      hsTot[i] = _hsTemp[i]; // + _hform[i];
      
      
      
      
//       CFout << "getSpeciesTotEnthalpies Ends" << "\n";
    }
	
	//CFout << "getspecies hsTot " << hsTot << "\n";
//     CFout << "getSpeciesTotEnthalpies temperature" << temp <<"\n";
//     CFout <<"getSpeciesTotEnthalpies -- > enthalpy: " << hTot << "\n";
	
// CFdouble massTot=0.0;
// 
//   for (CFint is = 0; is < _NS; ++is) {
//     cf_always_assert (_ys[is] > 1.1);
// 
//     const CFreal mm = _ys[is]/_mmasses[is];
//     massTot += mm;
//     _xs[is] = mm;
//   }
//   _xs *= 1./massTot;
	
    }
  
  }

  
  
  
//   Simple CALLING to EnthPoly
// /*/*  CFint i;//, imol = 0;
//   for(i = 0; i <_NS ;i++) {
//    
//     EnthPoly(temp,i);
//     
//     hsTot[i] = _hsP[i];*/*/
    
//     if(_flagMoleculesIDs[i])
//           {
//             (*hsVib)[imol]=0.0;
//             imol++;
//           }


    

    
//   }



//////////////////////////////////////////////////////////////////////////////

void CombustionModelLibrary::transportCoeffNEQ(CFreal& temp, //it was temperature, before
                                         CFdouble& pressure,
                                         CFreal* tVec,
                                         RealVector& normConcGradients,
                                         CFreal& eta,
                                         CFreal& lambdaTrRo,
                                         RealVector& lambdaInt,
                                         RealVector& rhoUdiff)
{
  lambdaTrRo=0.0;
  lambdaInt[0]=0.0;
  _massTtmTt = 0.0;
  CFdouble PHIs = 0.0;
  CFdouble lambda = 0.0;
  CFdouble CpP = 0.0;

  
//   CFdouble Cp = 0.0;
//   CFdouble rho = 0.0;
//   CFdouble eta = 0.0;
  
  
//   _D = 0.0;
  
//Viscosity function is not called, but rewritten here, because here calling eta() as function
// gives compiling error
//  CFout << "TransportCoeffNEQ Begins" << "\n";
  
  _EpsKb = 337.6241;
  _Tstar = 0.003*temp;
  _Ov = 1.16145*std::pow(_Tstar, -0.14874)+0.52487*exp(-0.77320*_Tstar)+2.16178*exp(-2.43787*_Tstar);
  _Fc = 1-0.2756*_AF;
  
  // C4H6 Butadiene requires specific viscosity model; it is placed as first member of species viscosity
  // array (place [0]) - _mmasses[0] is converted inside the calculation from kg to g again
  _NIUsTemp[0]=40.785*((_Fc*sqrt(_mmasses[0]*1000*temp))/(std::pow(_Vc,0.6667)*_Ov));
  _NIUsTemp[0]=_NIUsTemp[0]*std::pow(10.,-7.);
   // CFout << _NIUsTemp[0] << "\n";
  

  _massTtmS[0]=_ys[0]/_mmasses[0];
  
  
  switch(_chemID) 
  {
  //CFout << "Temperature in transportCoeffNEQ" << temp << "\n";  
  case 0:
  //CFout <<"Switching to eta() function case 0 \n";
  CFint i;//,j;
  for(i=1;i< _NS;i++)
    {
        _massTtmS[i]=_ys[i]/_mmasses[i];
	if (i==1){
	  if (temp <= 373.2) {
	     _NIUsTemp[i]=1.2277*std::pow(10.,-5.);
// 	     ofstream viscfile((_libPath+"validation/visc1.txt").c_str());
// 	     viscfile.open("visc1.txt");
// 	     CFout << _NIUsTemp[i] << "\n";
// 	     viscfile << "prova";
// 	     viscfile.close();
	     
	  }
	  else if (temp>373.2 && temp<1075.0) {
	          _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]);
		  _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
// 		   CFout << _NIUsTemp[i] << "\n";
	  }
	  else {
	       _NIUsTemp[i]=exp(_AsMu2[i]*log(temp)+_BsMu2[i]/temp+_CsMu2[i]/(temp*temp)+_DsMu2[i]);
	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
// 	        CFout << _NIUsTemp[i] << "\n";
	  } 
	}
	else if (i != 1) {
	   if (temp < 1000) {
	      _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]);
	      _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
// 	       CFout << _NIUsTemp[i] << "\n";
	  }
	  else {
	       _NIUsTemp[i]=exp(_AsMu2[i]*log(temp)+_BsMu2[i]/temp+_CsMu2[i]/(temp*temp)+_DsMu2[i]);

	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
// 	        CFout << _NIUsTemp[i] << "\n";
	  }
	}
      }

  
  break;
      
  case 1:
  //CFout <<"Switching to eta() function case 1 \n";
  //CFout << "Temperature in transportCoeffNEQ" << temp << "\n";
 
  for(i=1;i< _NS;i++)
    {
        _massTtmS[i]=_ys[i]/_mmasses[i];
	//H2O HAS to be the species in position [1] of the vectors (.cmix file accordingly!)
	if (i==1){
	  if (temp <= 373.2) {
	     _NIUsTemp[i]=1.2277*std::pow(10.,-5.);
	  }
	    else if (temp>373.2 && temp<=1075.0) {
	    _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]);
	    _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	  else {
	       _NIUsTemp[i]=exp(_AsMu2[i]*log(temp)+_BsMu2[i]/temp+_CsMu2[i]/(temp*temp)+_DsMu2[i]);
	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  } 
	}
	else if (i != 1) {
	   if (temp < 1000) {
	      _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]);
	      _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	  else {
	       _NIUsTemp[i]=exp(_AsMu2[i]*log(temp)+_BsMu2[i]/temp+_CsMu2[i]/(temp*temp)+_DsMu2[i]);
	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	}
      }
  // C4H6 Butadiene requires specific viscosity model; it is placed as first member of species viscosity
  // array (place [0])

  
  break;
  
  case 2:
  
  //  CFout << "Temperature in transportCoeffNEQ" << temp << "\n";
  //CFout <<"Switching to eta() function case 2 \n";
  for(i=1;i< _NS;i++)
    {
        _massTtmS[i]=_ys[i]/_mmasses[i];
	//H2O HAS to be the species in position [1] of the vectors (.cmix file accordingly!)
	if (i==1){
	  if (temp <= 373.2) {
	     _NIUsTemp[i]=1.2277*std::pow(10.,-5.);
	  }
	    else if (temp>373.2 && temp<=1075.0) {
	    _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]);
	    _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	  else {
	       _NIUsTemp[i]=exp(_AsMu2[i]*log(temp)+_BsMu2[i]/temp+_CsMu2[i]/(temp*temp)+_DsMu2[i]);
	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  } 
	}
	//H HAS to be the species in position [2] of the vectors (.cmix file accordingly!)
	else if (i==2) {
	  if (temp <=1000.0) {
	     _NIUsTemp[i]=1.4236*std::pow(10.,-5.);
	  }  
	  else {
	       _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]); 
	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	}
	  //O HAS to be the species in position [2] of the vectors (.cmix file accordingly!)
	else if (i==3) {
	  if (temp <=1000.0) {
	     _NIUsTemp[i]=4.9965*std::pow(10.,-5.);
	  }
	  else {
	       _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]); 
	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	}
	  // This can probably be substituted by a simple "else" without if
	else { //if (i != 1 && i != 2 && i != 3) {
	   if (temp < 1000) {
	      _NIUsTemp[i]=exp(_AsMu1[i]*log(temp)+_BsMu1[i]/temp+_CsMu1[i]/(temp*temp)+_DsMu1[i]); 
	      _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	
	  else {
	       _NIUsTemp[i]=exp(_AsMu2[i]*log(temp)+_BsMu2[i]/temp+_CsMu2[i]/(temp*temp)+_DsMu2[i]);
	       _NIUsTemp[i]=_NIUsTemp[i]*std::pow(10.,-7.);
	  }
	}
      }
  // C4H6 Butadiene requires specific viscosity model; it is placed as first member of species viscosity
  // array (place [0])

   break;
    
   }
   
 
CFint i,j;
    // Mixture viscosity calculation as sum weighted over PHI.
    for(i=0;i< _NS;i++)
    {
        PHIs=0.0;
	// PHIs is the same for ATDModelLibrary and CombustionModelLibrary
        for(j=0;j<_NS;j++)
        {
            PHIs+=_massTtmS[j]*std::pow(1+sqrt(_NIUsTemp[i]/_NIUsTemp[j])*std::pow(_mmasses[j]/_mmasses[i],0.25),2)/sqrt(8*(1+_mmasses[i]/_mmasses[j]));
        }

        eta+=_massTtmS[i]*_NIUsTemp[i]/PHIs;
    }  
  
  lambda = lambdaNEQ(temp, pressure);
  
  CpP = CpPoly(temp);
  
  for ( i = 0; i < _NS; ++i) {
     
	 rhoUdiff[i]=-normConcGradients[i]*(lambda)/(CpP*_Le);
	//It was: 
    //rhoUdiff[i]=-normConcGradients[i]*(lambda*_Rgas)/(CpP*_Le);
    
  }  
  

 CFout << "TransportCoeffNEQ Ends" << "\n";

}

//////////////////////////////////////////////////////////////////////////////


void CombustionModelLibrary::ReadDataMixture()
{
  string line = "start"; // Before was string line; not working on banchieri
  CFint linecount=0;

  _NS=0;

  ifstream mixfile((m_libPath+"mixtures/comb/"+_mixtureName+".cmix").c_str());

  if (mixfile.is_open())
  {

    while ( line != "stop" ) // Before was ( mixfile.good() ), not working on banchieri
    {
      getline(mixfile,line);

      // code to interpret the lines
      if (line[0] != '/')
      {
        linecount+=1;
        if (linecount==1)
        {
          _NS = static_cast<CFuint> (atoi(line.c_str()));
	  
	  //
	  //PREVIOUS
	  //
          
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

	  //
	  // NEW
	  //
    _xs.resize(_NS);  
    _AsMu1.resize(_NS);
    _BsMu1.resize(_NS);
    _CsMu1.resize(_NS);
    _DsMu1.resize(_NS);
    _AsMu2.resize(_NS);
    _BsMu2.resize(_NS);
    _CsMu2.resize(_NS);
    _DsMu2.resize(_NS);
    _AsK1.resize(_NS);
    _BsK1.resize(_NS);
    _CsK1.resize(_NS);
    _DsK1.resize(_NS);
    _AsK2.resize(_NS);
    _BsK2.resize(_NS);
    _CsK2.resize(_NS);
    _DsK2.resize(_NS);
    _a1cp1.resize(_NS);
    _a2cp1.resize(_NS);
    _a3cp1.resize(_NS);
    _a4cp1.resize(_NS);
    _a5cp1.resize(_NS);
    _a1cp2.resize(_NS);
    _a2cp2.resize(_NS);
    _a3cp2.resize(_NS);
    _a4cp2.resize(_NS);
    _a5cp2.resize(_NS);
    _a1h1.resize(_NS);
    _a2h1.resize(_NS);
    _a3h1.resize(_NS);
    _a4h1.resize(_NS);
    _a5h1.resize(_NS);
    _a6h1.resize(_NS);
    _a1h2.resize(_NS);
    _a2h2.resize(_NS);
    _a3h2.resize(_NS);
    _a4h2.resize(_NS);
    _a5h2.resize(_NS);
    _a6h2.resize(_NS);
    _KappasTemp.resize(_NS);
    _CpTemp.resize(_NS);
    _hsTemp.resize(_NS);
    _hsP.resize(_NS);
    _omegaTempDot.resize(_NS);
    
    
    
//
// CONTROL if _thermoID is needed!!!!!
//

         switch (_thermoID){
             case 0:
               _numTCR = 11;
               _numTC = 22;
               _thermoCoefs.resize(_numTC*_NS);
               break;
//             case 1:
//               _numTCR=6;
//               _numTC=30;
//               _thermoCoefs.resize(_numTC*_NS);
//               break;
//             case 2:
//               _numTCR=8;
//               _numTC=24;
//               _thermoCoefs.resize(_numTC*_NS);
//               break;
           }

// NUMBER INSIDE PARENTHESES HAS TO BE VERIFIED WITH CHEMISTRY FOLDERS
//           switch (_chemID){
//             case 0:
//               _chemCoefsTemp.resize(8);
//               break;
//             case 1:
//               _chemCoefsTemp.resize(8);
//               break;
//             case 2:
//               _chemCoefsTemp.resize(8);
//               break;
// 
//           }

        }
        if (linecount>1 && linecount<_NS+2)  _speciesNames[linecount-2]=line;
      }
    }
    mixfile.close();
  }

  else {
    CFout << "Unable to open mixture file: " << m_libPath+"mixtures/comb/"+_mixtureName+".cmix" << "\n";
    abort();
  }

}

//////////////////////////////////////////////////////////////////////////////
/*
void ATDModelLibrary::ReadDataChem()
{
  
}*/

//////////////////////////////////////////////////////////////////////////////

void CombustionModelLibrary::ReadDataSpecies(CFint i)
{
  string line = "start"; // Before was string line; not working on banchieri
  CFint linecount=0;
  CFdouble MMaux;
  stringstream StrStream;

// Three switches are created in order to read from three different folders, one for each chemical 
//  model, each one with different number of species' file to be read
  
  switch (_chemID){
    
   case 0:
   {
   // Changed the file extension to read from .spc (ATDModelLibrary) to .cspc (CombustionModelLibrary)
   CFout << "ReadDataSpecies case 0 - Venkateswaran. \n";
   ifstream spcfile((m_libPath+"species/venka/"+_speciesNames[i]+".cspc").c_str());
 
   if (spcfile.is_open())
   {
     while ( line != "stop" ) // Before was ( mixfile.good() ), not working on banchieri
	 {
       getline (spcfile,line);
       //CFout << "Species opened. \n";

       // code to interpret the lines
       if (line[0] != '/')
       { 
	 //CFout << "If Opened. \n";
         StrStream.str(line);
         linecount++;
         if ((linecount==1 && (StrStream>>_CA[i]).fail()) || (linecount==2 && (StrStream>>MMaux).fail())
             || (linecount==3 && (StrStream>>_hform[i]).fail()) || (linecount==4 && (StrStream>>_tau_vs[i]).fail())
             || (linecount==5 && (StrStream>>_Tets[i]).fail()) || (linecount==6 && (StrStream>>_AsMu1[i]).fail())
             || (linecount==7 && (StrStream>>_BsMu1[i]).fail()) || (linecount==8 && (StrStream>>_CsMu1[i]).fail())
	     || (linecount==9 && (StrStream>>_DsMu1[i]).fail()) || (linecount==10 && (StrStream>>_AsMu2[i]).fail())
	     || (linecount==11 && (StrStream>>_BsMu2[i]).fail()) || (linecount==12 && (StrStream>>_CsMu2[i]).fail())
	     || (linecount==13 && (StrStream>>_DsMu2[i]).fail()) || (linecount==14 && (StrStream>>_AsK1[i]).fail())
	     || (linecount==15 && (StrStream>>_BsK1[i]).fail()) || (linecount==16 && (StrStream>>_CsK1[i]).fail())
	     || (linecount==17 && (StrStream>>_DsK1[i]).fail()) || (linecount==18 && (StrStream>>_AsK2[i]).fail())
	     || (linecount==19 && (StrStream>>_BsK2[i]).fail()) || (linecount==20 && (StrStream>>_CsK2[i]).fail())
	     || (linecount==21 && (StrStream>>_DsK2[i]).fail()))
         {
           CFout << "Wrong format of species file - failure inside StrStream if: "<< m_libPath+"species/venka/"+_speciesNames[i]+".cspc" << "\n";
           abort();
         }

         StrStream.clear();
       }
     }
     _mmasses[i]=MMaux*1e-3;  // conversion from g to kg
     spcfile.close();
   }
   else {
     CFout << "Unable to open mixture file: " << m_libPath+"species/venka/"+_speciesNames[i]+".cspc" << "\n";
     abort();
  }

// IF I LEAVE THIS, IT STOPS - WHY???
//   if (linecount!=21)
//   {
//       CFout << "Wrong format of species file - failure in linecount if: " << m_libPath+"species/venka/"+_speciesNames[i]+".cspc" << "\n";
//       abort();
//   }

  break;
  }    
   case 1:
   {
   // Changed the file extension to read from .spc (ATDModelLibrary) to .cspc (CombustionModelLibrary)
   CFout << "ReadDataSpecies case 1 - Jones-Lindstedt 4 reactions. \n";
   ifstream spcfile((m_libPath+"species/jones4/"+_speciesNames[i]+".cspc").c_str());
 
   if (spcfile.is_open())
   {
     while ( line != "stop" ) // Before was ( mixfile.good() ), not working on banchieri
	 {
       getline (spcfile,line);

       // code to interpret the lines
       if (line[0] != '/')
       {
         StrStream.str(line);
         linecount++;
         if ((linecount==1 && (StrStream>>_CA[i]).fail()) || (linecount==2 && (StrStream>>MMaux).fail())
             || (linecount==3 && (StrStream>>_hform[i]).fail()) || (linecount==4 && (StrStream>>_tau_vs[i]).fail())
             || (linecount==5 && (StrStream>>_Tets[i]).fail()) || (linecount==6 && (StrStream>>_AsMu1[i]).fail())
             || (linecount==7 && (StrStream>>_BsMu1[i]).fail()) || (linecount==8 && (StrStream>>_CsMu1[i]).fail())
	     || (linecount==9 && (StrStream>>_DsMu1[i]).fail()) || (linecount==10 && (StrStream>>_AsMu2[i]).fail())
	     || (linecount==11 && (StrStream>>_BsMu2[i]).fail()) || (linecount==12 && (StrStream>>_CsMu2[i]).fail())
	     || (linecount==13 && (StrStream>>_DsMu2[i]).fail()) || (linecount==14 && (StrStream>>_AsK1[i]).fail())
	     || (linecount==15 && (StrStream>>_BsK1[i]).fail()) || (linecount==16 && (StrStream>>_CsK1[i]).fail())
	     || (linecount==17 && (StrStream>>_DsK1[i]).fail()) || (linecount==18 && (StrStream>>_AsK2[i]).fail())
	     || (linecount==19 && (StrStream>>_BsK2[i]).fail()) || (linecount==20 && (StrStream>>_CsK2[i]).fail())
	     || (linecount==21 && (StrStream>>_DsK2[i]).fail()))
         {
           CFout << "Wrong format of species file - failure inside StrStream if: "<< m_libPath+"species/jones4/"+_speciesNames[i]+".cspc" << "\n";
           abort();
         }

         StrStream.clear();
       }
     }
     _mmasses[i]=MMaux*1e-3;  // conversion from g to kg
     spcfile.close();
   }
   else {
     CFout << "Unable to open mixture file: " << m_libPath+"species/jones4/"+_speciesNames[i]+".cspc" << "\n";
     abort();
  }

// IF I LEAVE THIS, IT STOPS - WHY???
//   if (linecount!=21)
//   {
//       CFout << "Wrong format of species file - failure in linecount if: " << m_libPath+"species/jones4/"+_speciesNames[i]+".cspc" << "\n";
//       abort();
//   }

  break;
  }
   case 2:
   {
   // Changed the file extension to read from .spc (ATDModelLibrary) to .cspc (CombustionModelLibrary)
   CFout << "ReadDataSpecies case 2 - Jones-Lindstedt 6 reactions. \n";
   ifstream spcfile((m_libPath+"species/jones6/"+_speciesNames[i]+".cspc").c_str());
 
   if (spcfile.is_open())
   {
     while ( line != "stop" ) // Before was ( mixfile.good() ), not working on banchieri
	 {
       getline (spcfile,line);

       // code to interpret the lines
       if (line[0] != '/')
       {
         StrStream.str(line);
         linecount++;
         if ((linecount==1 && (StrStream>>_CA[i]).fail()) || (linecount==2 && (StrStream>>MMaux).fail())
             || (linecount==3 && (StrStream>>_hform[i]).fail()) || (linecount==4 && (StrStream>>_tau_vs[i]).fail())
             || (linecount==5 && (StrStream>>_Tets[i]).fail()) || (linecount==6 && (StrStream>>_AsMu1[i]).fail())
             || (linecount==7 && (StrStream>>_BsMu1[i]).fail()) || (linecount==8 && (StrStream>>_CsMu1[i]).fail())
	     || (linecount==9 && (StrStream>>_DsMu1[i]).fail()) || (linecount==10 && (StrStream>>_AsMu2[i]).fail())
	     || (linecount==11 && (StrStream>>_BsMu2[i]).fail()) || (linecount==12 && (StrStream>>_CsMu2[i]).fail())
	     || (linecount==13 && (StrStream>>_DsMu2[i]).fail()) || (linecount==14 && (StrStream>>_AsK1[i]).fail())
	     || (linecount==15 && (StrStream>>_BsK1[i]).fail()) || (linecount==16 && (StrStream>>_CsK1[i]).fail())
	     || (linecount==17 && (StrStream>>_DsK1[i]).fail()) || (linecount==18 && (StrStream>>_AsK2[i]).fail())
	     || (linecount==19 && (StrStream>>_BsK2[i]).fail()) || (linecount==20 && (StrStream>>_CsK2[i]).fail())
	     || (linecount==21 && (StrStream>>_DsK2[i]).fail()))
         {
           CFout << "Wrong format of species file - failure inside StrStream if: "<< m_libPath+"species/jones6/"+_speciesNames[i]+".cspc" << "\n";
           abort();
         }

         StrStream.clear();
       }
     }
     _mmasses[i]=MMaux*1e-3;  // conversion from g to kg
     spcfile.close();
   }
   else {
     CFout << "Unable to open mixture file: " << m_libPath+"species/jones6/"+_speciesNames[i]+".cspc" << "\n";
     abort();
  }

// IF I LEAVE THIS, IT STOPS - WHY???
//   if (linecount!=21)
//   {
//       CFout << "Wrong format of species file - failure in linecount if: " << m_libPath+"species/jones6/"+_speciesNames[i]+".cspc" << "\n";
//       abort();
//   }

  break;
 }
    
 }

}

//////////////////////////////////////////////////////////////////////////////

void CombustionModelLibrary::ReadDataThermo(CFint i)
{
  string line = "start"; // Before was string line; not working on banchieri
  int linecount=0;
  stringstream StrStream;

// Changed the file extension to read from .thm (ATDModelLibrary) to .cthm (CombustionModelLibrary)
  
  ifstream thmfile((m_libPath+"thermo/comb/"+_speciesNames[i]+_thermoName+".cthm").c_str());

  if (thmfile.is_open())
  {
    //CFint iind=i*_numTC;

    while ( line != "stop" ) // Before was ( mixfile.good() ), not working on banchieri
		{
      getline (thmfile,line);

      // code to interpret the lines

if (line[0] != '/')
             {
         StrStream.str(line);
         linecount++;
         if ((linecount==1 && (StrStream>>_a1cp1[i]).fail()) || (linecount==2 && (StrStream>>_a2cp1[i]).fail())
             || (linecount==3 && (StrStream>>_a3cp1[i]).fail()) || (linecount==4 && (StrStream>>_a4cp1[i]).fail())
             || (linecount==5 && (StrStream>>_a5cp1[i]).fail()) || (linecount==6 && (StrStream>>_a1cp2[i]).fail())
             || (linecount==7 && (StrStream>>_a2cp2[i]).fail()) || (linecount==8 && (StrStream>>_a3cp2[i]).fail())
	     || (linecount==9 && (StrStream>>_a4cp2[i]).fail()) || (linecount==10 && (StrStream>>_a5cp2[i]).fail())
	     || (linecount==11 && (StrStream>>_a1h1[i]).fail()) || (linecount==12 && (StrStream>>_a2h1[i]).fail())
	     || (linecount==13 && (StrStream>>_a3h1[i]).fail()) || (linecount==14 && (StrStream>>_a4h1[i]).fail())
	     || (linecount==15 && (StrStream>>_a5h1[i]).fail()) || (linecount==16 && (StrStream>>_a6h1[i]).fail())
	     || (linecount==17 && (StrStream>>_a1h2[i]).fail()) || (linecount==18 && (StrStream>>_a2h2[i]).fail())
	     || (linecount==19 && (StrStream>>_a3h2[i]).fail()) || (linecount==20 && (StrStream>>_a4h2[i]).fail())
	     || (linecount==21 && (StrStream>>_a5h2[i]).fail()) || (linecount==22 && (StrStream>>_a6h2[i]).fail()))
         {
           CFout << "Wrong format of thermo file:"<< m_libPath+"thermo/comb/"+_speciesNames[i]+_thermoName+".cthm" << "\n";
          abort();
         }

         StrStream.clear();
       }
     }
// Old version of thermofile reading
// //       if (line[0] != '/')
// //       {
// //         StrStream.str(line);
// //         linecount++;
// //         if (linecount>_numTC)
// //         {
// //             CFout << "Too many coefficients on file:"<< m_libPath+"thermo/comb/"+_speciesNames[i]+_thermoName+".cthm" << "\n";
// //             abort();
// //         }
// //         if ((StrStream>>_thermoCoefs[iind+linecount-1]).fail())
// //         {
// //           CFout << "Wrong format of thermo file:"<< m_libPath+"thermo/comb/"+_speciesNames[i]+_thermoName+".cthm" << "\n";
// //           abort();
// //         }
// //         StrStream.clear();
// //       }
// //     }
    thmfile.close();
  }
  else {
    CFout << "Unable to open mixture file" << "thermo/comb/"+_speciesNames[i]+_thermoName+".cthm" << "\n";
    abort();
  }

  if (linecount<_numTC)
  {
      CFout << "Not enough coefficients on file:"<< m_libPath+"thermo/comb/"+_speciesNames[i]+_thermoName+".cthm" << "\n";
      abort();
  }
}

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

void CombustionModelLibrary:: setLibrary()
{
  
CFint i;//,j;
// constant values 
_AF=0.192;
_Tc=425.17;
_Vc=220;
_EpsKb=337.6241;

  //universal constants
  _NA=6.02214129e23; //Avogadro constant from NIST - Codata 2010
  _KB=1.3806488e-23; //Boltzmann constant from NIST - Codata 2010
  _PI=3.14159265;
  _Rgas=_NA*_KB;
  _Tref=298.16;
  _Le=1.0; //Le = 1, first guess

//   //At the present stage we are only considering the two-temperature model and neutral species
//   _nbTvib=1;
//   _nbTe=0;

  CFout << "CombustionModelLibrary::libpath     = "  << m_libPath << "\n";
  CFout << "CombustionModelLibrary::mixtureName = "  << _mixtureName << "\n";
  CFout << "CombustionModelLibrary::chemName = " << _chemName << "\n";
//   CFout << "CombustionModelLibrary::thermoName = " << _thermoName << "\n";

  // if nobody has configured the library path, it is done here
  // with default = basedir + "plugins/ATDModel/data/"
  if (m_libPath == "") {
    CFLog(NOTICE, "CombustionModelLibrary::libpath set to default" << "\n");
    std::string baseDir = Environment::DirPaths::getInstance().getBaseDir().string();
    m_libPath = baseDir + "/plugins/ATDModel/data/";
  }

// Define _thermoID
   if (_thermoName=="polynomial") _thermoID=0;
//   _thermoID definition remains in order to permit future implementation of further
//   thermal models (e.g. "direct" spline)
//   else if (_thermoName=="Gnoffo") _thermoID=1;
//   else if (_thermoName=="McBride") _thermoID=2;
   else {
       CFout << "ThermoName" + _thermoName + "is not valid." << "\n";
       abort();
   }


// Three chemical systems are defined
// Define _chemID

if (_chemName=="venka2r") _chemID=0;
else if (_chemName=="jones4r") _chemID=1;
else if (_chemName=="jones6r") _chemID=2;
else {
     CFout << "chemName " + _chemName + " is not valid." << "\n";
     abort();
     
}


// //   Define _chemID
//   if (_chemName=="air5park1" || _chemName=="nitrogen2park1") _chemID=0;
//   else if (_chemName=="air5park2" || _chemName=="nitrogen2park2") _chemID=1;
//   else if (_chemName=="air5park3" || _chemName=="nitrogen2park3") _chemID=2;
//   else if (_chemName=="air5DunnKang" || _chemName=="nitrogen2DunnKang") _chemID=3;
//   else {
//       CFout << "chemName " + _chemName + " is not valid." << "\n";
//       abort();
//   }

//   // Test _TempID
//   if (_TempID!=0 && _TempID!=1) {
//       CFout << "TempID " << _TempID << " is not valid." << "\n";
//       abort();
//   }

//   // Test _RelaxID
//   if (_RelaxID!=0 && _RelaxID!=1 && _RelaxID!=2 && _RelaxID!=3) {
//       CFout << "RelaxID " << _RelaxID << " is not valid." << "\n";
//       abort();
//   }

//   // Test _TetsID
//   if (_TetsID!=0 && _TetsID!=1) {
//       CFout << "TetsID " << _TetsID << " is not valid." << "\n";
//       abort();
//   }

  //Read data from files initialize arrays
  ReadDataMixture();

//   ReadDataChem();

CFout << "Mixture data read. \n";

  //Read data from files
  for (i = 0; i < _NS; ++i) {
     
     ReadDataThermo(i);
     CFout << "Thermo data read. \n";
 
     ReadDataSpecies(i);
     CFout << "Species data read. \n";
     
  }

//   setMoleculesIDs(_molecIDs);

//   _flagMoleculesIDs.resize(_NS,false);
//   for (CFuint i = 0; i < _molecIDs.size(); ++i) {
//     _flagMoleculesIDs[_molecIDs[i]] = true;
//   }
//   for (CFint i = 0; i < _NS; ++i) {
//     _atomicityCoeff[i] = (_flagMoleculesIDs[i]) ? 2.5 : 1.5;
//   }

//   for(i=0;i<_NS;i++)
//   {
//       for(j=0;j<_NS;j++)
//       {
//           _NIUsr(i,j)=1000*(_mmasses[i]*_mmasses[j])/(_mmasses[i]+_mmasses[j]);
//           _Asr(i,j)=(_tau_vs[i]!=-1)? 1.16e-3*sqrt(_NIUsr(i,j))*std::pow(_tau_vs[i],4.0/3):0;
//       }
//   }

  CFout << "Library set." << "\n";
  
  
  
// //VALIDATION FOR THERMOPHYSICAL AND TRANSPORT PROPERTIES
// //SPECIES VISCOSITY VALIDATION
// // viscfile and condfile has to be rearranged every time in order to be consistent with the number of species 
// // (AND THEIR ORDER!!!) in the following validation print-to-file dummy calls.
// // ALSO [i] value is to be changed accordingly to species order inside chemical models  
  
// // ofstream viscfile1("validation/viscC4H6.txt");
// // ofstream viscfile2("validation/viscH2O.txt");
// // ofstream viscfile3("validation/viscH.txt");
// // ofstream viscfile4("validation/viscO.txt");
// // ofstream viscfile5("validation/viscCO.txt");
// // ofstream viscfile6("validation/viscCO2.txt");
// // ofstream viscfile7("validation/viscOH.txt");
// // ofstream viscfile8("validation/viscO2.txt");
// // ofstream viscfile9("validation/viscH2.txt");
// //  
//   CFdouble temp;
//   CFdouble lambda;
//   CFdouble pressure = 1.0;
//   for(temp = 300; temp<=5000; temp++)
//   {
// 
//   lambda=lambdaNEQ(temp, pressure);
//   
//   CFout << "lambda: " << lambda <<"\n";
//   CFout << "_KappasTemp: " << _KappasTemp <<"\n";
//   
//   }
// // 
// //   ofstream viscfile1("validation/viscC4H6.txt",ios::app);
// // 
// // if (viscfile1.is_open()) 
// // {
// //   viscfile1 << _NIUsTemp[0] <<"\n";
// //   viscfile1.close();
// // }
// // 
// //   ofstream viscfile2("validation/viscH2O.txt",ios::app);
// // 
// // if (viscfile2.is_open()) 
// // {
// //   viscfile2 << _NIUsTemp[1] <<"\n";
// //   viscfile2.close();
// // }
// // 
// //   ofstream viscfile3("validation/viscH.txt",ios::app);
// // 
// // if (viscfile3.is_open()) 
// // {
// //   viscfile3 << _NIUsTemp[2] <<"\n";
// //   viscfile3.close();
// // }
// // 
// //   ofstream viscfile4("validation/viscO.txt",ios::app);
// // 
// // if (viscfile4.is_open()) 
// // {
// //   viscfile4 << _NIUsTemp[3] <<"\n";
// //   viscfile4.close();
// // }
// // 
// //   ofstream viscfile5("validation/viscCO.txt",ios::app);
// // 
// // if (viscfile5.is_open()) 
// // {
// //   viscfile5 << _NIUsTemp[4] <<"\n";
// //   viscfile5.close();
// // }
// // 
// //   ofstream viscfile6("validation/viscCO2.txt",ios::app);
// // 
// // if (viscfile6.is_open()) 
// // {
// //   viscfile6 << _NIUsTemp[5] <<"\n";
// //   viscfile6.close();
// // }
// // 
// //   ofstream viscfile7("validation/viscOH.txt",ios::app);
// // 
// // if (viscfile7.is_open()) 
// // {
// //   viscfile7 << _NIUsTemp[6] <<"\n";
// //   viscfile7.close();
// // }
// // 
// //   ofstream viscfile8("validation/viscO2.txt",ios::app);
// // 
// // if (viscfile8.is_open()) 
// // {
// //   viscfile8 << _NIUsTemp[7] <<"\n";
// //   viscfile8.close();
// // }
// // 
// //   ofstream viscfile9("validation/viscH2.txt",ios::app);
// // 
// // if (viscfile9.is_open()) 
// // {
// //   viscfile9 << _NIUsTemp[8] <<"\n";
// //   viscfile9.close();
// // }
// // 
// //   }
// //   
// // }
  
// ofstream hfile1("validation/hC4H6.txt");
// ofstream hfile2("validation/hH2O.txt");
// ofstream hfile3("validation/hO2.txt");
// //ofstream hfile4("validation/hH2.txt");
// ofstream hfile4("validation/hCO.txt");
// ofstream hfile5("validation/hCO2.txt");
// ofstream hfile7("validation/hOH.txt");
// ofstream hfile8("validation/hO2.txt");
// ofstream hfile9("validation/hH2.txt");
// ofstream viscfile1("validation/viscC4H6prova.txt");
// ofstream condfile2("validation/condH2O.txt");
// ofstream condfile3("validation/condO2.txt");
// ofstream condfile4("validation/condH2.txt");
// ofstream condfile4("validation/condCO.txt");
// ofstream condfile5("validation/condCO2.txt");
// ofstream condfile7("validation/condOH.txt");
// ofstream condfile8("validation/condO2.txt");
// ofstream condfile9("validation/condH2.txt");
 
//   CFdouble temp;
//   for(temp = 300; temp<=5000; temp++)
//   {
//    EnthPoly(temp);
// // 
//    ofstream hfile1("validation/hC4H6.txt",ios::app);
// 
//   if (hfile1.is_open()) 
//    {
//     hfile1 << _hsTemp[0] <<"\n";
//     hfile1.close();
//    }
//    
//    ofstream hfile2("validation/hH2O.txt",ios::app);
// 
//   if (hfile2.is_open()) 
//    {
//     hfile2 << _hsTemp[1] <<"\n";
//     hfile2.close();
//    }
// 
//    ofstream hfile3("validation/hO2.txt",ios::app);
// 
//   if (hfile3.is_open()) 
//    {
//     hfile3 << _hsTemp[2] <<"\n";
//     hfile3.close();
//    }

//    ofstream hfile4("validation/hH2.txt",ios::app);
// 
//   if (hfile4.is_open()) 
//    {
//     hfile4 << _hsTemp[3] <<"\n";
//     hfile4.close();
//    }

//    ofstream hfile4("validation/hCO.txt",ios::app);
// 
//   if (hfile4.is_open()) 
//    {
//     hfile4 << _hsTemp[3] <<"\n";
//     hfile4.close();
//    }
// 
// 
//    ofstream hfile5("validation/hCO2.txt",ios::app);
// 
//   if (hfile5.is_open()) 
//    {
//     hfile5 << _hsTemp[4] <<"\n";
//     hfile5.close();
//    }
//    
//       ofstream hfile7("validation/hOH.txt",ios::app);
// 
//   if (hfile7.is_open()) 
//    {
//     hfile7 << _hsTemp[6] <<"\n";
//     hfile7.close();
//    }
//    
//       ofstream hfile8("validation/hO2.txt",ios::app);
// 
//   if (hfile8.is_open()) 
//    {
//     hfile8 << _hsTemp[7] <<"\n";
//     hfile8.close();
//    }
//    
//       ofstream hfile9("validation/hH2.txt",ios::app);
// 
//   if (hfile9.is_open()) 
//    {
//     hfile9 << _hsTemp[8] <<"\n";
//     hfile9.close();
//    }


//   }
//   
//   ofstream condfile2("validation/condH2O.txt",ios::app);
//   
//    if (condfile2.is_open()) 
//   {
//    condfile2 << _KappasTemp[1] <<"\n";
//    condfile2.close();
//   }
//    
//    ofstream condfile3("validation/condO2.txt",ios::app);
//    
//    if (condfile3.is_open()) 
//   {
//    condfile3 << _KappasTemp[2] <<"\n";
//    condfile3.close();
//   }
//   
// //    ofstream condfile4("validation/condH2.txt",ios::app);
// //   
// //    if (condfile4.is_open()) 
// //   {
// //    condfile4 << _KappasTemp[3] <<"\n";
// //    condfile4.close();
// //   }
//   
//   
//   ofstream condfile4("validation/condCO.txt",ios::app);
//   
//    if (condfile4.is_open()) 
//   {
//    condfile4 << _KappasTemp[3] <<"\n";
//    condfile4.close();
//   }
//   
//   ofstream condfile5("validation/condCO2.txt",ios::app);
//   
//    if (condfile5.is_open()) 
//   {
//    condfile5 << _KappasTemp[4] <<"\n";
//    condfile5.close();
//   }
//   
// //   ofstream viscfile1("validation/viscC4H6prova.txt",ios::app);
// //  if (viscfile1.is_open()) 
// //   {
// //    viscfile1 << _NIUsTemp[0] <<"\n";
// //    viscfile1.close();
// //   }
//   
//   }

  

//   ofstream condfile2("validation/condH2O.txt",ios::app);
// 
// if (condfile2.is_open()) 
// {
//   condfile2 << _KappasTemp[1] <<"\n";
//   condfile2.close();
// }
// 
//   ofstream condfile3("validation/condH.txt",ios::app);
// 
// if (condfile3.is_open()) 
// {
//   condfile3 << _KappasTemp[2] <<"\n";
//   condfile3.close();
// }
// 
//   ofstream condfile4("validation/condO.txt",ios::app);
// 
// if (condfile4.is_open()) 
// {
//   condfile4 << _KappasTemp[3] <<"\n";
//   condfile4.close();
// }
// 
//   ofstream condfile5("validation/condCO.txt",ios::app);
// 
// if (condfile5.is_open()) 
// {
//   condfile5 << _KappasTemp[4] <<"\n";
//   condfile5.close();
// }
// 
//   ofstream condfile6("validation/condCO2.txt",ios::app);
// 
// if (condfile6.is_open()) 
// {
//   condfile6 << _KappasTemp[5] <<"\n";
//   condfile6.close();
// }
// 
//   ofstream condfile7("validation/condOH.txt",ios::app);
// 
// if (condfile7.is_open()) 
// {
//   condfile7 << _KappasTemp[6] <<"\n";
//   condfile7.close();
// }
// 
//   ofstream condfile8("validation/condO2.txt",ios::app);
// 
// if (condfile8.is_open()) 
// {
//   condfile8 << _KappasTemp[7] <<"\n";
//   condfile8.close();
// }
// 
//   ofstream condfile9("validation/condH2.txt",ios::app);
// 
// if (condfile9.is_open()) 
// {
//   condfile9 << _KappasTemp[8] <<"\n";
//   condfile9.close();
// }

  
//   ofstream cpfile1("validation/cpC4H6.txt");
// 
//   CFdouble temp;
//   for(temp = 300; temp<=5000; temp++)
//   {
// 
//   if (temp <= 1000) {
// 	     _CpTemp[0]= _a1cp1[0]*std::pow(temp,4)+_a2cp1[0]*std::pow(temp,3)+_a3cp1[0]*std::pow(temp,2)+_a4cp1[0]*temp+_a5cp1[0];
// 	  }
// 	  else {
// 	       _CpTemp[0]= _a1cp2[0]*std::pow(temp,4)+_a2cp2[0]*std::pow(temp,3)+_a3cp2[0]*std::pow(temp,2)+_a4cp2[0]*temp+_a5cp2[0];
// 	  } 
// 
//   ofstream cpfile1("validation/cpC4H6.txt",ios::app);
// 
//  if (cpfile1.is_open()) 
//   {
//    cpfile1 << _CpTemp[0] <<"\n";
//    cpfile1.close();
//   }
//   
//   CFdouble temp, rho;
// //const RealVector ys;
//   RealVector omega;
// 
// // ATDModelLibrary::setSpeciesFractions();
// 
// RealVector ys;
// CFout << ys << "\n";
// //   
//   for(temp = 300; temp < 5000;temp++){
//    
//     CFout << "inside temp loop" << "\n";
//    
//     for(rho=0.1;rho<2;rho = rho + 0.1) {
//       
//       CFout << "inside rho loop" << "\n";
//   
// 	for(i=1;i<_NS;i++) {
// 	
// 	  CFout << "inside i loop" << "\n";
// 	  
// 	  omegaContribution(temp,rho,_ys,omega);
// 	  
// 	}
//     
//     }
//   
//   }


// }

}


    //

//////////////////////////////////////////////////////////////////////////////

} // namespace ATDModel

} // namespace Physics

} // namespace COOLFluiDs

//////////////////////////////////////////////////////////////////////////////

