#include "MutationUsage.hh"
#include "ComputeTransportCoefs.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/Stopwatch.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::Mutation;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MutationUsage {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ComputeTransportCoefs, ComputeWithMutation, MutationUsageModule, 1>
computeTransportCoefsProv("TrCoefs");

//////////////////////////////////////////////////////////////////////////////

void ComputeTransportCoefs::compute()
{
 cout << "### ComputeTransportCoefs::compute() ###" << endl;

 SafePtr<MutationLibrary> library = this->getLibrary();

 double p = 101325.0;
 cout << "enter pressure" << endl;
 cin >> p;
 cout << "pressure = " << p << endl;

 const double Tmin = 10.0; // !!!
 const double Tmax = 10000.0;
 const double deltaT = 10.0;

 // Random elemental fraction Ye
 const CFuint nbElements = library->getNbElements();
 const CFuint nbSpecies = library->getNbSpecies();
 RealVector Ye(nbElements);
 RealVector Xe(nbElements);

 CFdouble lambdaR; // Thermal reactive conductivity
 CFdouble lambdaD; // Thermal demixing conductivitu
 RealVector lambdaEL(nbElements); // Elemental heat transfer coefficients
 RealMatrix eldifcoef(nbElements, nbElements); // Elemental multicomponent diffusion coefficients
 RealVector eltdifcoef(nbElements); // Elemental thermal demixing coefficients

// CFreal sumRand = 0.;
//
// for(CFint i = 0; i < nbElements; ++i) {
//   Ye[i] = (CFdouble(rand())/RAND_MAX);
//   sumRand += Ye[i];
// }
// for(CFint i = 0; i < nbElements; ++i) {
//   Ye[i] /= sumRand;
//   cout << "Ye[" << i << "] = " << Ye[i] << " ";
// }
// cout << endl;
 Ye[0] = 0.0;
 Ye[1] = 0.0;

 Stopwatch<> stp;
 stp.start();

 // Set elemental fractions if you do not want to use Mutation defaults
// library->setElemFractions(Ye);

 const std::string outfile = "mutation_TrCoefs_p" +
   StringOps::to_str(static_cast<CFuint>(p)) + ".dat";

 ofstream fout(outfile.c_str());
 RealVector Xs(nbSpecies);
 RealVector Ys(nbSpecies);
 RealVector Omega(nbSpecies);
//  bool flg_jac = false;
 RealMatrix Jacobian(nbSpecies+2,nbSpecies+2);
 double T = Tmin;


 // Writing Tecplot header
 fout << "VARIABLES=T ";
 for(CFuint i = 0; i < nbElements; ++i) {
  fout << ",Ye" << i << " ";
 }
 for(CFuint i = 0; i < nbSpecies; ++i) {
  fout << ",Xs" << i << " ";
 }
 for(CFuint i = 0; i < nbSpecies; ++i) {
  fout << ",Ys" << i << " ";
 }
 fout << ",lambdR ,lambdD ";
 for(CFuint i = 0; i < nbElements; ++i) {
  fout << ",lambdEL" << i << " ";
 }
 fout << ",lambdELTot "; // For gases with 2 elements
 for(CFuint i = 0; i < nbElements; ++i) {
  for(CFuint j = 0; j < nbElements; ++j) {
   fout << ",DTot" << i << j << " ";
  }
 }
 for(CFuint i = 0; i < nbElements; ++i) {
  fout << ",DT" << i << " ";
 }
  for(CFuint i = 0; i < nbSpecies; ++i) {
  fout << ",Omega" << i << " ";
 }
 fout << endl;

 do {
//    CFdouble rho = library->density(T,p,CFNULL);
   library->setComposition(T,p,&Xs);
   library->getTransportCoefs(T, p,
			      lambdaR, lambdaD, lambdaEL,
			      eldifcoef, eltdifcoef);
//   CFdouble TNEQ = T+400; // + 2000.0; // Can move out of equilibrium
// TO BE FIXED MARCO
//   RealVector TVIB2[0] = 400; // + 2000.0; // Can move out of equilibrium
//   library->getMassProductionTerm(TNEQ, TVIB2, p, rho, Ys, flg_jac, Omega, Jacobian);

   fout << T << " ";
   for(CFuint i = 0; i < nbElements; ++i) {
    fout << Ye[i] << " ";
   }
   for(CFuint i = 0; i < nbSpecies; ++i) {
    fout << Xs[i] << " ";
   }
   for(CFuint i = 0; i < nbSpecies; ++i) {
    fout << Ys[i] << " ";
   }
   fout << lambdaR << " " << lambdaD << " ";
   for(CFuint i = 0; i < nbElements; ++i) {
    fout << lambdaEL[i] << " ";
   }
   fout << lambdaEL[1]-lambdaEL[0] << " "; // For gases with 2 elements
   for(CFuint i = 0; i < nbElements; ++i) {
    for(CFuint j = 0; j < nbElements; ++j) {
     fout << 2.*eldifcoef(i,j) << " "; // For gases with 2 elements (2.*)
    }
   }
   for(CFuint i = 0; i < nbElements; ++i) {
    fout << eltdifcoef[i] << " ";
   }
   for(CFuint i = 0; i < nbSpecies; ++i) {
    fout << Omega[i] << " ";
   }
   fout << endl;

   T += deltaT;
 } while (T < Tmax);

 cout << "time: " << stp.read() << " sec" << endl;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MutationUsage

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
