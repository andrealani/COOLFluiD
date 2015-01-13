#include "GreyRadiator.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalModel.hh"

namespace COOLFluiD {

namespace RadiativeTransfer {

Environment::ObjectProvider<GreyRadiator,
                            Radiator,
                            RadiativeTransferModule,
                            1>
GreyRadiatorProvider("GreyRadiator");

Environment::ObjectProvider<GreyWallRadiator,
                            Radiator,
                            RadiativeTransferModule,
                            1>
GreyWallRadiatorProvider("GreyWallRadiator");


void GreyRadiator::defineConfigOptions(Config::OptionList& options)
{
  options.addConfigOption< CFreal >("ElemEmsCoeff",
          "Element Emission coeff (W /sr / m^(2,3) /Angstrom /m)");
  options.addConfigOption< CFreal >("ElemAbsCoeff",
          "Element Absorption coeff (CellKappa for a volume; WallAlpha for a surface)");
  options.addConfigOption< bool >("allIsGrey",
          "if all computational domain spectral properties are constant");
  options.addConfigOption< CFuint >("nbComPlankTerms",
          "Number of terms of the series to be calculated for the Comulative Plank Function");
}

GreyRadiator::GreyRadiator(const std::string& name):
     Radiator(name)
{
    addConfigOptionsTo(this);

    m_emsCoeff=-1.;
    setParameter("ElemEmsCoeff",&m_emsCoeff);

    m_absCoeff=-1.;
    setParameter("ElemAbsCoeff",&m_absCoeff);

    m_allIsGrey = false;
    setParameter("allIsGrey",&m_allIsGrey);

    m_nbComPlankTerms = 10;
    setParameter("nbComPlankTerms",&m_nbComPlankTerms);

}

void GreyRadiator::configure(Config::ConfigArgs& args){
    ConfigObject::configure(args);
}

//////////////////////////////////////////////////////////////////////////////

CFreal GreyRadiator::computeComulativePlankFraction(const CFreal lambda,
                                                    const CFreal T      ){
  //calculate the partial Sum
  const CFreal x=m_c2/(lambda*T);
  CFreal sum=0.;
  for(CFreal n=1.; n<CFreal(m_nbComPlankTerms); ++n){
    sum+=(x*x*x/n + 3.*x*x/(n*n) + 6.*x/(n*n*n)+ 6./(n*n*n*n) )*exp(-n*x);
  }
 // std::cout<<" lambda: "<<lambda<<" m_c1: "<< m_c1 <<"m_nbComPlankTerms: "<<m_nbComPlankTerms<<" sum: "<<sum<<std::endl;
  return m_c1*sum;
}

//////////////////////////////////////////////////////////////////////////////

inline CFreal GreyRadiator::computePlank(const CFreal lambda, const CFreal T){
  return m_c3/(pow(lambda,5) * ( std::exp(m_c2/(lambda*T) ) -1) );
}

inline CFreal GreyRadiator::computeStefanBoltzmann(const CFreal T){
  return m_StefanBoltzmann*pow(T,4);
}

//////////////////////////////////////////////////////////////////////////////

inline CFreal GreyRadiator::getCurrentCellTemperature(){
    CFuint stateID = m_radPhysicsHandlerPtr->getCurrentCellStateID();
    static Framework::DataHandle<Framework::State*, Framework::GLOBAL> m_states
            = m_radPhysicsHandlerPtr->getDataSockets()->states.getDataHandle();

    //std::cout<<"temperature: "<<(&(*m_states[stateID])[0])[m_tempID]<<std::endl;
    return (&(*m_states[stateID])[0])[m_tempID];
}

inline CFreal GreyRadiator::getCurrentWallTemperature(){
    CFuint wallGeoIdx = m_radPhysicsHandlerPtr->getCurrentWallTrsIdx();
    Framework::State* state = m_radPhysicsPtr->getWallState( wallGeoIdx );

//    CFuint wallGeoID = m_radPhysicsHandlerPtr->getCurrentWallGeoID();
//    Framework::DataHandle<CFreal> faceNormals
//      = m_radPhysicsHandlerPtr->getDataSockets()->normals.getDataHandle();
 //   std::cout<<"normal "<<faceNormals[DIM_2D*wallGeoID]<<' '<< faceNormals[DIM_2D*wallGeoID+1]<<std::endl;
    
    CFreal Temperature = (&(*state)[0])[m_tempID];
    cf_assert(Temperature > 0 && "Temperature negative, is the temp ID well defined?");
    return Temperature;
}

//////////////////////////////////////////////////////////////////////////////

CFreal GreyRadiator::getAbsorption( CFreal lambda, RealVector &s_o ){
  return m_absCoeff;
}

//////////////////////////////////////////////////////////////////////////////

CFreal GreyRadiator::getEmission( CFreal lambda, RealVector &s_o ){
  CFreal T = getCurrentElemTemperature();
  return m_emsCoeff * computePlank(lambda*m_angstrom, T);
}

//////////////////////////////////////////////////////////////////////////////

void GreyRadiator::setup(){
  m_tempID = m_radPhysicsHandlerPtr->getTempID();
  m_TRStypeID = m_radPhysicsPtr->getTRStypeID();
}

void GreyRadiator::setupSpectra(CFreal wavMin, CFreal wavMax)
{
 m_minWav = wavMin*m_angstrom;
 m_maxWav = wavMax*m_angstrom;
}

CFreal GreyRadiator::getSpectaLoopPower(){
   CFreal T = getCurrentElemTemperature();

   CFreal a1 = computeComulativePlankFraction(m_minWav, T);
   CFreal a2 = computeComulativePlankFraction(m_maxWav, T);
   CFreal a3 = (a2-a1)*computeStefanBoltzmann(T) * m_emsCoeff * getSpaceIntegrator();

//   if(m_TRStypeID == WALL ){
//      std::cout<<"Temperatute: "<<T<<" Area: "<<getSpaceIntegrator()<<" emsCoeff "<<m_emsCoeff<<std::endl;
//      std::cout<<"computeStefanBoltzmann "<<computeStefanBoltzmann(T)<<std::endl;
//      std::cout<<"space integrator:  "<<getSpaceIntegrator()<<std::endl;
//      std::cout<<"total:  "<<(a2-a1)*a3 <<std::endl;
//      std::cout<<"a1: "<<a1<<" a2: "<<a2<<" a3: "<<a3<<std::endl;
//   }

   return (a2-a1)*a3;
}

void GreyRadiator::getRandomEmission(CFreal &lambda, RealVector &s_o){
    CFreal T = getCurrentElemTemperature();

    static CFuint dim = Framework::PhysicalModelStack::getActive()->getDim();
    static CFuint dim2 = m_radPhysicsHandlerPtr->isAxi() ? 3 : dim;
    //m_TRStypeID=MEDIUM;
    //emission direction: diffuse
    //std::cout<<"dim1: "<<dim2<<" dim2_: "<<s_o.size()<<std::endl;
    getRandomDirections(dim2, s_o);

    if(!m_allIsGrey){
      //emission wavelength: bissection method
      //@TODO Maybe a better algorithm?
      bool test;
      CFreal a = m_minWav;
      CFreal b = m_maxWav;
      CFreal minComulativePlank =computeComulativePlankFraction(a,T);
      CFreal maxComulativePlank =computeComulativePlankFraction(b,T);

      CFreal c = (m_minWav+m_maxWav)/2.;
      CFreal target = m_rand.uniformRand();
      CFreal m_tol=1e-10;
      CFreal f_c;
      do{
        f_c =(computeComulativePlankFraction(c,T)-minComulativePlank)
            / (maxComulativePlank-minComulativePlank);
        test = f_c > target;
        b=test?c:b;
        a=test?a:c;
        c= (a+b)/2.;
      }  while ( abs(f_c - target) > m_tol );

      lambda = f_c * m_angstrom;
    }
    else{
      lambda = 1.;
    }
}

}
}
