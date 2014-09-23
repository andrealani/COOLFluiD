#include "MutationUsage.hh"
#include "ComputeRhoHTOverAt.hh"
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

Environment::ObjectProvider<ComputeRhoHTOverAt, ComputeWithMutation,  MutationUsageModule, 1>
computeRhoHTOverAtProv("RhoHTOverAt");

//////////////////////////////////////////////////////////////////////////////

void ComputeRhoHTOverAt::compute()
{
 cout << "### ComputeRhoHTOverAt::compute() ###" << endl;

 SafePtr<MutationLibrary> library = this->getLibrary();

 double p = 101325.0;
 cout << "enter pressure" << endl;
 cin >> p;
 cout << "pressure = " << p << endl;

 const double Tmin = 100.0;
 const double Tmax = 30000.0;
 const double deltaT = 50.0;

 Stopwatch<> stp;
 stp.start();

 const std::string outfile = "mutation_rhoHtOverAt_p" +
   StringOps::to_str(static_cast<CFuint>(p)) + ".dat";

 CFreal tol = 1e-6;
 ofstream fout(outfile.c_str());
 RealVector x(library->getNbSpecies());
 RealVector dhe(3);
 RealVector pertDhe(3);
 double T = Tmin;
 do {
   library->setComposition(T,p,&x);
   library->setDensityEnthalpyEnergy(T,p,dhe);

   CFdouble dT = tol*T;
   CFdouble pertT = T + dT;
   library->setComposition(pertT,p,&x);
   library->setDensityEnthalpyEnergy(pertT,p,pertDhe);

   const CFdouble rho = dhe[0];
   const CFdouble h = dhe[1];
   const CFdouble e = dhe[2];
   const CFdouble rhoT = (pertDhe[0] - dhe[0])/dT;
   const CFdouble hT = (pertDhe[1] - dhe[1])/dT;
   const CFdouble eT = (pertDhe[2] - dhe[2])/dT;
   cout << "rhoT = " << rhoT << endl;
   cout << "hT = " << hT << endl;
   cout << "eT = " << eT << endl;

   const CFdouble coeff = rho*hT/(rhoT*e - rhoT*h + rho*eT);

   fout << T << " " << coeff << endl;
   T += deltaT;
 } while (T < Tmax);

 cout << "time: " << stp.read() << " sec" << endl;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MutationUsage

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
