#include "FiniteVolumeMHD/FiniteVolumeMHD.hh"

#include "FiniteVolumeMHD/MHD3DInitStatePFSS.hh"
#include "Common/CFLog.hh"
#include "Framework/MethodCommandProvider.hh"
#include "Framework/DataHandle.hh"
#include "Common/BadValueException.hh"
#include "Common/FilesystemException.hh"
#include "MHD/MHD3DVarSet.hh"

#include <gsl/gsl_sf_legendre.h>

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Physics::MHD;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Numerics {

    namespace FiniteVolume {

//////////////////////////////////////////////////////////////////////////////

MethodCommandProvider<MHD3DInitStatePFSS, CellCenterFVMData, FiniteVolumeMHDModule>
mHD3DInitStatePFSSProvider("MHD3DInitStatePFSS");

//////////////////////////////////////////////////////////////////////////////

void MHD3DInitStatePFSS::defineConfigOptions(Config::OptionList& options)
{
   options.addConfigOption< CFreal >("mass","Mass of the external object for which the escape velocity will be computed.");
   options.addConfigOption< CFreal >("TRef","Reference temperature to compute the isothermal speed of sound for Parker's solar wind.");
   options.addConfigOption< CFreal >("epsilon","Desired accuracy for Parker solution.");
   options.addConfigOption< CFreal >("rMin","Radius of the inner boundary.");
   options.addConfigOption< CFreal >("rMax","Radius of the outer boundary.");
   options.addConfigOption< CFreal >("rSource","Radius of the source surface for the PFSS model.");
   options.addConfigOption< CFreal >("n","Polytropic index.");
}

//////////////////////////////////////////////////////////////////////////////

MHD3DInitStatePFSS::MHD3DInitStatePFSS(const std::string& name) :
  InitState(name),
  _varSet(CFNULL),
  socket_Almreal("Almreal"),
  socket_Almimg("Almimg"),
  socket_Blmreal("Blmreal"),
  socket_Blmimg("Blmimg"),
  _namePFSSDataFile() 
{
  addConfigOptionsTo(this);

  // the mass of the Sun is assigned by default
  _mass = 1.98892e30;
  setParameter("mass",&_mass);

  _TRef = 1.5e6;
  setParameter("TRef",&_TRef);

  _epsilon = 1.0e-6;
  setParameter("epsilon",&_epsilon);

  _rMin = 1.03;
  setParameter("rMin",&_rMin);

  _rMax = 20.0;
  setParameter("rMax",&_rMax);

  _rSource = 2.5;
  setParameter("rSource",&_rSource);

  _n = 1.05;
  setParameter("n",&_n);

}

//////////////////////////////////////////////////////////////////////////////

MHD3DInitStatePFSS::~MHD3DInitStatePFSS()
{
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DInitStatePFSS::readInputFile()
{
  DataHandle<std::vector<CFreal> > Almreal  = socket_Almreal.getDataHandle();
  DataHandle<std::vector<CFreal> > Almimg  = socket_Almimg.getDataHandle();
  DataHandle<std::vector<CFreal> > Blmreal  = socket_Blmreal.getDataHandle();
  DataHandle<std::vector<CFreal> > Blmimg  = socket_Blmimg.getDataHandle();

  // The order of coefficients should be Almreal, Almimg, Blmreal and Blmimg each separated with a blank line in the input file 

  SelfRegistPtr<Environment::FileHandlerInput> fhandle = Environment::SingleBehaviorFactory<Environment::FileHandlerInput>::getInstance().create();
  ifstream& inputFile = fhandle->open(constructFilename(_namePFSSDataFile));

  std::string line, tmpstr;
  vector<CFreal> tmprow;

  const CFuint nbLModes = _varSet->getNbLModes();

  for (CFuint iLine = 0; iLine <= nbLModes; ++iLine) {
      getline(inputFile,line);
      istringstream liness( line );
  
      for (CFuint jLine = 0; jLine <= nbLModes; ++jLine) {
       	getline(liness, tmpstr, ' ');
       	tmprow.push_back(atof(tmpstr.c_str())); 
      }
      Almreal[iLine] = tmprow;
      tmprow.clear();
  }

  // skip the undesired modes and the blank line separating the coefficients
  while (getline(inputFile,line)) {
        if (line == "") {
           break;
        }
  }

  for (CFuint iLine = 0; iLine <= nbLModes; ++iLine) {
      getline(inputFile,line);
      istringstream liness( line );

      for (CFuint jLine = 0; jLine <= nbLModes; ++jLine) {
        getline(liness, tmpstr, ' ');
        tmprow.push_back(atof(tmpstr.c_str()));
      }
      Almimg[iLine] = tmprow;
      tmprow.clear();
  }

  // skip the undesired modes and the blank line separating the coefficients
  while (getline(inputFile,line)) {
        if (line == "") {
           break;
        }
  }

  for (CFuint iLine = 0; iLine <= nbLModes; ++iLine) {
      getline(inputFile,line);
      istringstream liness( line );

      for (CFuint jLine = 0; jLine <= nbLModes; ++jLine) {
        getline(liness, tmpstr, ' ');
        tmprow.push_back(atof(tmpstr.c_str()));
      }
      Blmreal[iLine] = tmprow;
      tmprow.clear();
  }

  // skip the undesired modes and the blank line separating the coefficients
  while (getline(inputFile,line)) {
        if (line == "") {
           break;
        }
  }

  for (CFuint iLine = 0; iLine <= nbLModes; ++iLine) {
      getline(inputFile,line);
      istringstream liness( line );

      for (CFuint jLine = 0; jLine <= nbLModes; ++jLine) {
        getline(liness, tmpstr, ' ');
        tmprow.push_back(atof(tmpstr.c_str()));
      }
      Blmimg[iLine] = tmprow;
      tmprow.clear();
  }

  fhandle->close();
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DInitStatePFSS::computePFSSMagneticField(const RealVector& stateCoordsSpherical,
					RealVector& BPFSSCartesian,
                                	RealMatrix& sphCarTransMat)
{
  CFreal r = stateCoordsSpherical[0];
  const CFreal theta = stateCoordsSpherical[1];
  const CFreal phi = stateCoordsSpherical[2];
  const CFreal rTemp = r;

  if (r >= _rSource) {
     r = _rSource;
  }

  DataHandle<std::vector<CFreal> > Almreal  = socket_Almreal.getDataHandle();
  DataHandle<std::vector<CFreal> > Almimg  = socket_Almimg.getDataHandle();
  DataHandle<std::vector<CFreal> > Blmreal  = socket_Blmreal.getDataHandle();
  DataHandle<std::vector<CFreal> > Blmimg  = socket_Blmimg.getDataHandle();

  const CFuint nbLModes = _varSet->getNbLModes();

  CFreal Br = 0.0, Btheta = 0.0, Bphi = 0.0;

  for (CFuint l = 0; l <= nbLModes; ++l) {
      for (CFuint m = 0; m <= l; ++m) {
          // Assign the real part of B(r,phi)          
          const CFreal Ylmampl = gsl_sf_legendre_sphPlm(l,m,cos(theta));
          const CFreal Brreal = -(Almreal[l])[m]*(double)l*pow(r,(double)l-1.0)+(Blmreal[l])[m]*((double)l+1.0)*pow(r,-(double)l-2.0);
          const CFreal Brimg = -(Almimg[l])[m]*(double)l*pow(r,(double)l-1.0)+(Blmimg[l])[m]*((double)l+1.0)*pow(r,-(double)l-2.0);
          const CFreal BthetaCoeff = -1.0/(r*sin(theta));
          const CFreal Bphireal = BthetaCoeff*(double)m*((Almreal[l])[m]*pow(r,(double)l)+(Blmreal[l])[m]*pow(r,-(double)l-1.0));
          const CFreal Bphiimg = BthetaCoeff*(double)m*((Almimg[l])[m]*pow(r,(double)l)+(Blmimg[l])[m]*pow(r,-(double)l-1.0));

          const CFreal Brampl = sqrt(Brreal*Brreal+Brimg*Brimg);
          const CFreal Bphiampl = sqrt(Bphireal*Bphireal+Bphiimg*Bphiimg);

	  // TBD: In order to obtain the correct polarity of the magnetic dipoles on the photosphere -= is used instead of += for Br, Bphi and Btheta. This should be checked.
          const CFreal Brangle = atan2(Brimg,Brreal);
	  Br -= Brampl*Ylmampl*cos(Brangle+(double)m*phi);

          const CFreal Bphiangle = atan2(Bphiimg,Bphireal);
	  Bphi -= Bphiampl*Ylmampl*cos(Bphiangle+(double)m*phi+0.5*MathTools::MathConsts::CFrealPi());
      }
  }

  for (CFuint l = 1; l <= (nbLModes-1); ++l) {
      for (CFuint m = 0; m <= l; ++m) {
          const CFreal Rlm = sqrt(((double)l*(double)l-(double)m*(double)m)/(4.0*(double)l*(double)l-1.0));
          const CFreal Rlpl1m = sqrt((((double)l+1.0)*((double)l+1.0)-(double)m*(double)m)/(4.0*((double)l+1.0)*((double)l+1.0)-1.0));
          // Assign the real part of B(theta)          
          const CFreal Ylmampl = gsl_sf_legendre_sphPlm(l,m,cos(theta));
          const CFreal BthetaCoeff = -1.0/(r*sin(theta));
          const CFreal Bthetareal = BthetaCoeff*(Rlm*((double)l-1.0)*((Almreal[l-1])[m]*pow(r,(double)l-1.0)+(Blmreal[l-1])[m]*pow(r,-(double)l))
                                            -Rlpl1m*((double)l+2.0)*((Almreal[l+1])[m]*pow(r,(double)l+1.0)+(Blmreal[l+1])[m]*pow(r,-(double)l-2.0)));
          const CFreal Bthetaimg = BthetaCoeff*(Rlm*((double)l-1.0)*((Almimg[l-1])[m]*pow(r,(double)l-1.0)+(Blmimg[l-1])[m]*pow(r,-(double)l))
                                            -Rlpl1m*((double)l+2.0)*((Almimg[l+1])[m]*pow(r,(double)l+1.0)+(Blmimg[l+1])[m]*pow(r,-(double)l-2.0)));

          const CFreal Bthetaampl = sqrt(Bthetareal*Bthetareal+Bthetaimg*Bthetaimg);

          const CFreal Bthetaangle = atan2(Bthetaimg,Bthetareal);
	  Btheta -= Bthetaampl*Ylmampl*cos(Bthetaangle+(double)m*phi);
      }
  }

  RealVector BPFSSSpherical(PhysicalModelStack::getActive()->getDim());

  BPFSSSpherical[0] = Br;
  if (rTemp < _rSource) {
        BPFSSSpherical[1] = Btheta;
        BPFSSSpherical[2] = Bphi;
  } else {
        if (rTemp > _rSource) {
              // B field is radial and decreases by r^2 beyond the source surface
              BPFSSSpherical[0] = Br/(rTemp*rTemp);
        }
        BPFSSSpherical[1] = 0.0;
        BPFSSSpherical[2] = 0.0;
  }

  BPFSSCartesian = sphCarTransMat*BPFSSSpherical;
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DInitStatePFSS::computeParkerSolution(const CFreal r, 
					  RealVector& velParkerSpherical,
					  CFreal& rhoParker,
					  CFreal& pParker)
{
  // gravitational constant
  const CFreal G = 6.67384e-11;

  // Boltzmann constant
  const CFreal k = 1.3806503e-23;

  // magnetic permeability at vacuum
  const CFreal mu0 = 4.e-7*MathTools::MathConsts::CFrealPi();

  // mass of proton
  const CFreal mp = 1.67262158e-27;

  // mass of electron
  const CFreal me = 9.10938188e-31;

  const CFreal nRef = _varSet->getNRef();
  const CFreal BRef = _varSet->getBRef();
  const CFreal lRef = _varSet->getLRef();

  const CFreal rhoRef = nRef*(mp+me);
  const CFreal vRef = BRef / sqrt(mu0*rhoRef);

  // calculate the escape velocity: sqrt(2GM/R)
  const CFreal vEsc = sqrt(2.0*G*_mass/lRef)/vRef;

  // calculate the isothermal speed of sound
  const CFreal a = sqrt(2.0*k*_TRef/mp)/vRef;

  // calculate the expected transsonic point
  const CFreal rTrans = (vEsc*vEsc)/(4.0*a*a);

  // check whether the transsonic point is in the computational domain
  cf_assert(rTrans > _rMin);
  cf_assert(rTrans < _rMax);
  cf_assert(r != rTrans);

  // Parker's solar wind is radial
  velParkerSpherical[1] = 0.0;
  velParkerSpherical[2] = 0.0; 

  CFreal ur0 = 1.0, ur1, urh, diff;
  bool hasSolution = false;

  // calculate the Parker's solar wind speed that goes from subsonic to supersonic regime
  if (r < rTrans) {
	for (CFuint iter = 1; iter < 1000; ++iter) {	
	    //ur1 = ((vEsc*vEsc*vEsc*vEsc)/(16.0*r*r))*exp(0.5*((ur0*ur0)+3.0-((vEsc*vEsc)/r)));
	    ur1 = ((rTrans*rTrans)/(r*r))*exp(0.5*((ur0*ur0)+3.0-(4.0*rTrans/r)));
            diff = fabs(ur1-ur0);
	    if (diff < _epsilon) {
		velParkerSpherical[0] = ur1*a;
                hasSolution = true;
                break;
            }
	    else {
	        ur0 = ur1;
	    }
     	}
        cf_assert(hasSolution);
  }
  if (r > rTrans) {
        for (CFuint iter = 1; iter < 1000; ++iter) {
            //urh = ((vEsc*vEsc)/r)-3.0+2.0*log(16.0*ur0*r*r/(vEsc*vEsc*vEsc*vEsc));
	    urh = (4.0*rTrans/r)-3.0+2.0*log(ur0*r*r/(rTrans*rTrans));
	    ur1 = sqrt(urh);
            diff = fabs(ur1-ur0);
            if (diff < _epsilon) {
                velParkerSpherical[0] = ur1*a;
                hasSolution = true;
                break;
            }
            else {
                ur0 = ur1;
            }
        }
        cf_assert(hasSolution);
  }
  
  // calculate the approximate base solar wind speed in the vicinity of the photosphere (i.e. at _rMin)
  const CFreal vBase = a*((rTrans*rTrans)/(_rMin*_rMin))*exp(1.5-(2.0*rTrans/_rMin));
  // density obeys rho*v_r*r^2 = const.
  // assuming rhoBase = 1
  rhoParker = (vBase*_rMin*_rMin)/(velParkerSpherical[0]*r*r);

  // relation between density and pressure is p=a*a*rho
  // assuming rhoBase = 1
  const CFreal pBase = a*a;

  // a polytropic relation is assumed between density and pressure, p/rho^n = const.
  pParker = pBase*pow(rhoParker,_n);
}

//////////////////////////////////////////////////////////////////////////////

void MHD3DInitStatePFSS::executeOnTrs()
{
  SafePtr<TopologicalRegionSet> trs = getCurrentTRS();
  CFLogDebugMax( "InitState::executeOnTrs() called for TRS: "
  << trs->getName() << "\n");

  if (trs->getName() != "InnerFaces") {
    throw BadValueException (FromHere(),"InitState not applied to InnerFaces!!!");
  }

  // this cannot be used for FV boundary faces because
  // ghost state and inner state could have the same local ID !!!
  SafePtr<vector<CFuint> > trsStates = trs->getStatesInTrs();
  DataHandle < Framework::State*, Framework::GLOBAL > states = socket_states.getDataHandle();

  std::vector<CFuint>::iterator itd;
  CFreal rhoParker, pParker;
  RealVector stateCoordsCartesian(PhysicalModelStack::getActive()->getDim()),
        stateCoordsSpherical(PhysicalModelStack::getActive()->getDim()),
        BPFSSCartesian(PhysicalModelStack::getActive()->getDim()),
        velParkerSpherical(PhysicalModelStack::getActive()->getDim()),
        velParkerCartesian(PhysicalModelStack::getActive()->getDim());
  RealMatrix carSphTransMat(PhysicalModelStack::getActive()->getDim(),PhysicalModelStack::getActive()->getDim()),
        sphCarTransMat(PhysicalModelStack::getActive()->getDim(),PhysicalModelStack::getActive()->getDim());
  const CFreal gamma = _varSet->getModel()->getGamma();
  /*const CFreal lRef = _varSet->getModel()->getLRef();
  const CFreal BRef = _varSet->getModel()->getBRef();*/
  if(_inputAdimensionalValues)
  {
    for (itd = trsStates->begin(); itd != trsStates->end(); ++itd) {
      State* const currState = states[(*itd)];
      stateCoordsCartesian = currState->getCoordinates();
      //stateCoordsCartesian *= lRef; 
      _varSet->setTransformationMatrices(stateCoordsCartesian,stateCoordsSpherical,carSphTransMat,sphCarTransMat);
      computePFSSMagneticField(stateCoordsSpherical,BPFSSCartesian,sphCarTransMat);
      computeParkerSolution(stateCoordsSpherical[0]/*/lRef*/,velParkerSpherical,rhoParker,pParker);
      velParkerCartesian = sphCarTransMat*velParkerSpherical;
      //_vFunction.evaluate(currState->getCoordinates(), *_input);
      // input and update variables are assumed to be conservative
      (*_input)[0] = rhoParker;
      (*_input)[1] = rhoParker*velParkerCartesian[0];
      (*_input)[2] = rhoParker*velParkerCartesian[1];
      (*_input)[3] = rhoParker*velParkerCartesian[2];
      (*_input)[4] = BPFSSCartesian[0]/*/BRef*/;
      (*_input)[5] = BPFSSCartesian[1]/*/BRef*/;
      (*_input)[6] = BPFSSCartesian[2]/*/BRef*/;
      (*_input)[7] = pParker/(gamma-1.0)+(0.5*((*_input)[1]*(*_input)[1]+(*_input)[2]*(*_input)[2]+
		(*_input)[3]*(*_input)[3])/(*_input)[0])+0.5*((*_input)[4]*(*_input)[4]+(*_input)[5]*(*_input)[5]+
		(*_input)[6]*(*_input)[6]);
      *currState = *_inputToUpdateVar->transform(_input);
    }
  }
  else
  {
    State dimState;
    for (itd = trsStates->begin(); itd != trsStates->end(); ++itd) {
      State* const currState = states[(*itd)];
      stateCoordsCartesian = currState->getCoordinates();
      //stateCoordsCartesian *= lRef; 
      _varSet->setTransformationMatrices(stateCoordsCartesian,stateCoordsSpherical,carSphTransMat,sphCarTransMat);
      computePFSSMagneticField(stateCoordsSpherical,BPFSSCartesian,sphCarTransMat);
      computeParkerSolution(stateCoordsSpherical[0]/*/lRef*/,velParkerSpherical,rhoParker,pParker);
      velParkerCartesian = sphCarTransMat*velParkerSpherical;
      //_vFunction.evaluate(currState->getCoordinates(), *_input);
      // input and update variables are assumed to be conservative
      (*_input)[0] = rhoParker;
      (*_input)[1] = rhoParker*velParkerCartesian[0];
      (*_input)[2] = rhoParker*velParkerCartesian[1];
      (*_input)[3] = rhoParker*velParkerCartesian[2];
      (*_input)[4] = BPFSSCartesian[0]/*/BRef*/;
      (*_input)[5] = BPFSSCartesian[1]/*/BRef*/;
      (*_input)[6] = BPFSSCartesian[2]/*/BRef*/;
      (*_input)[7] = pParker/(gamma-1.0)+(0.5*((*_input)[1]*(*_input)[1]+(*_input)[2]*(*_input)[2]+
                (*_input)[3]*(*_input)[3])/(*_input)[0])+0.5*((*_input)[4]*(*_input)[4]+(*_input)[5]*(*_input)[5]+
                (*_input)[6]*(*_input)[6]);
      dimState = *_inputToUpdateVar->transform(_input);
      _varSet->setAdimensionalValues(dimState, *currState);
    }
  }
}


//////////////////////////////////////////////////////////////////////////////

void MHD3DInitStatePFSS::setup()
{

  InitState::setup();

  _varSet = getMethodData().getUpdateVar().d_castTo<MHD3DVarSet>();
  cf_assert(_varSet.isNotNull());

  DataHandle<std::vector<CFreal> > Almreal  = socket_Almreal.getDataHandle();
  DataHandle<std::vector<CFreal> > Almimg  = socket_Almimg.getDataHandle();
  Almreal.resize(_varSet->getNbLModes()+1);
  Almimg.resize(_varSet->getNbLModes()+1);

  DataHandle<std::vector<CFreal> > Blmreal  = socket_Blmreal.getDataHandle();
  DataHandle<std::vector<CFreal> > Blmimg  = socket_Blmimg.getDataHandle();
  Blmreal.resize(_varSet->getNbLModes()+1);
  Blmimg.resize(_varSet->getNbLModes()+1);

  _namePFSSDataFile = _varSet->getNamePFSSCoeffFile();

  // if this is a parallel simulation, only ONE process at a time reads the file
  if (PE::GetPE().IsParallel()) {
    
    PE::GetPE().setBarrier();

    for (CFuint i = 0; i < PE::GetPE().GetProcessorCount(); ++i) {
      if (i == PE::GetPE().GetRank ()) {
        readInputFile();
      }
      
      PE::GetPE().setBarrier();
    }
  }
  else {
    readInputFile();
  }
}

//////////////////////////////////////////////////////////////////////////////

boost::filesystem::path MHD3DInitStatePFSS::constructFilename(std::string fileName)
{
  boost::filesystem::path fpath(fileName);

  CFout << "Reading PFSS spherical harmonics coefficients from: " << fpath.string() << "\n";

  return fpath;
}

//////////////////////////////////////////////////////////////////////////////

std::vector<Common::SafePtr<BaseDataSocketSource> >
MHD3DInitStatePFSS::providesSockets()
{
  std::vector<Common::SafePtr<BaseDataSocketSource> > result =
          InitState::providesSockets();
  result.push_back(&socket_Almreal);
  result.push_back(&socket_Almimg);
  result.push_back(&socket_Blmreal);
  result.push_back(&socket_Blmimg);
  return result;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace FiniteVolume

  } // namespace Numerics

} // namespace COOLFluiD
