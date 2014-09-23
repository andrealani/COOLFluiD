#include "NEQ/NEQ.hh"
#include "Euler2DNEQRoeVinokurToConsInRef.hh"
#include "NavierStokes/EulerPhysicalModel.hh"
#include "Framework/PhysicalModel.hh"
#include "Environment/ObjectProvider.hh"
#include "Framework/PhysicalChemicalLibrary.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::MathTools;
using namespace COOLFluiD::Physics::NavierStokes;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace Physics {

    namespace NEQ {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<Euler2DNEQRoeVinokurToConsInRef, VarSetMatrixTransformer,
NEQModule, 1> euler2DNEQRoeVinokurToConsInRefProvider("Euler2DNEQRoeVinokurToConsInRef");

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQRoeVinokurToConsInRef::Euler2DNEQRoeVinokurToConsInRef
(Common::SafePtr<Framework::PhysicalModelImpl> model) :
  VarSetMatrixTransformer(model),
  _model(model->getConvectiveTerm().d_castTo<NEQTerm>()),
  _RiGas()
{
}

//////////////////////////////////////////////////////////////////////////////

Euler2DNEQRoeVinokurToConsInRef::~Euler2DNEQRoeVinokurToConsInRef()
{
}

//////////////////////////////////////////////////////////////////////////////

void Euler2DNEQRoeVinokurToConsInRef::setMatrixFromRef()
{
//   std::cout << "void Euler2DNEQRoeVinokurToConsInRef::setMatrixFromRef()" << std::endl;
  cf_assert(_model.isNotNull());
  const RealVector& linearData = _model->getPhysicalData();
//   std::cout << "linearData = " << linearData <<  std::endl;
  // find a way of storing this pointer (sdd setup() function)
  static Common::SafePtr<PhysicalChemicalLibrary> library =
    PhysicalModelStack::getActive()->getImplementor()->
    getPhysicalPropertyLibrary<PhysicalChemicalLibrary>();
  
  Common::SafePtr<PhysicalChemicalLibrary::ExtraData> eData = library->getExtraData();
  
  const CFuint nbSpecies = _model->getNbScalarVars(0);
  const CFuint firstSpecies = _model->getFirstScalarVar(0);
  const CFuint firstTv = _model->getFirstScalarVar(1);

  const CFuint uID  = nbSpecies;
  const CFuint vID  = nbSpecies+1;
  const CFuint eID  = nbSpecies+2;
  const CFuint evID = nbSpecies+3;
  
  //!!!!!!!!! To be done in a setup method:
  _alpha_bar.resize(nbSpecies);

  for (CFuint iSpecies = 0; iSpecies < nbSpecies; iSpecies++){
    
    _alpha_bar[iSpecies] =  eData->dP_Bar[iSpecies];
  }

   
  _beta_bar = eData->dP_Bar[nbSpecies];

//   std::cout << "_alpha_bar = " << _alpha_bar << "; _beta_bar =" << _beta_bar << std::endl;

  const CFreal onePlusBeta_bar = 1. + _beta_bar;
  const CFreal oneOverOnePlusBeta_bar = 1. / onePlusBeta_bar;
  const CFreal beta_barOverOnePlusBeta_bar = _beta_bar / onePlusBeta_bar;
  //!!!!!!!!!
  
  const CFreal rho = linearData[EulerTerm::RHO];
  const CFreal sqrtRho = sqrt(rho);
  const CFreal p = linearData[EulerTerm::P];
  const CFreal T = linearData[EulerTerm::T];
  const CFreal H = linearData[EulerTerm::H];

  const CFreal ev = linearData[firstTv];
  const CFreal u = linearData[EulerTerm::VX];
  const CFreal v = linearData[EulerTerm::VY];
  const CFreal V2 = linearData[EulerTerm::V]*linearData[EulerTerm::V];


  const CFreal sqrtRhoU = sqrtRho * u;
  const CFreal sqrtRhoV = sqrtRho * v;
  const CFreal sqrtRhoH = sqrtRho * H;
  const CFreal sqrtRhoEv = sqrtRho * ev;
  
//   const CFreal cvTr = eData->dEdT; // check this !!!!
//   const CFuint firstTv = _model->getFirstScalarVar(1);
//   const CFuint firstSpecies = _model->getFirstScalarVar(0);
//   const CFreal ev = linearData[firstTv];
//   const CFreal dd = rho*T*cvTr + p;
//   const CFreal Rb = p/(sqrtRho*T);
//   const CFreal c = sqrtRho*cvTr + Rb;
//   const CFreal c2 = c*c;
//   const CFreal dRhoEdzev = sqrtRho*p/dd;

  _RiGas.resize(nbSpecies);
  library->setRiGas(_RiGas);

  CFreal sum_YsAlphas_bar = 0.;
  for (CFuint iSpecies = 0; iSpecies < nbSpecies; ++iSpecies) {

    const CFreal  Ys = linearData[firstSpecies + iSpecies];
    sum_YsAlphas_bar +=  Ys * _alpha_bar[iSpecies];
  }
//   cout << "Euler2DNEQRoeVinokurToConsInRef::_alpha_bar = " << _alpha_bar << endl;
//   cout << "Euler2DNEQRoeVinokurToConsInRef::_beta_bar = " << _beta_bar << endl;
  
  for (CFuint iSpecies = 0; iSpecies < nbSpecies; ++iSpecies) {
        
    for (CFuint jSpecies = 0; jSpecies < nbSpecies; ++jSpecies) {

      const CFreal  Ys = linearData[firstSpecies + iSpecies];
      _transMatrix(iSpecies,jSpecies) = sqrtRho * Ys;

      if (iSpecies == jSpecies){
        _transMatrix(iSpecies,jSpecies) += sqrtRho;
      }
      
    }
    
    _transMatrix(uID,iSpecies) = sqrtRhoU;
    _transMatrix(vID,iSpecies) = sqrtRhoV;

    _transMatrix(eID,iSpecies) = oneOverOnePlusBeta_bar * sqrtRhoH + beta_barOverOnePlusBeta_bar * sqrtRhoEv - oneOverOnePlusBeta_bar * sqrtRho * (_alpha_bar[iSpecies] + sum_YsAlphas_bar);//sqrtRho/c2*(cH + cEv + cV);
    _transMatrix(evID,iSpecies) = sqrtRhoEv;
  }
  
  _transMatrix(uID,uID) = sqrtRho;
  _transMatrix(vID,vID) = sqrtRho;
  
  _transMatrix(eID,uID)  = beta_barOverOnePlusBeta_bar * sqrtRho * u;
  _transMatrix(eID,vID)  = beta_barOverOnePlusBeta_bar * sqrtRho * v;
  _transMatrix(eID,eID)  = oneOverOnePlusBeta_bar * sqrtRho;
  _transMatrix(eID,evID) = beta_barOverOnePlusBeta_bar * sqrtRho;
  
  _transMatrix(evID,evID) = sqrtRho;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace NEQ

  } // namespace Physics

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
