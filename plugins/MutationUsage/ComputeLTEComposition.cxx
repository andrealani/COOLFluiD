#include "ComputeLTEComposition.hh"
#include "Environment/ObjectProvider.hh"
#include "Common/Stopwatch.hh"
#include "MutationUsage.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Framework;
using namespace COOLFluiD::Physics::Mutation;
using namespace COOLFluiD::Common;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace MutationUsage {

//////////////////////////////////////////////////////////////////////////////

Environment::ObjectProvider<ComputeLTEComposition, ComputeWithMutation, MutationUsageModule, 1>
computeLTECompositionProv("xLTE");

//////////////////////////////////////////////////////////////////////////////

void ComputeLTEComposition::compute()
{
 cout << "### ComputeLTEComposition::compute() ###" << endl;

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

 const std::string outfile = "mutation_X_p" +
   StringOps::to_str(static_cast<CFuint>(p)) + ".dat";

 ofstream fout(outfile.c_str());
 RealVector x(library->getNbSpecies());
 double T = Tmin;

 // Writing Tecplot header
 fout << "VARIABLES=T ";
 const CFuint nbSpecies = library->getNbSpecies();
 for(CFuint i = 0; i < nbSpecies; ++i) {
  fout << ",X[" << i << "] ";
 }
 fout << endl;

 do {
   library->setComposition(T,p,&x);
   fout << T << " " << x << endl;
   T += deltaT;
 } while (T < Tmax);

 cout << "time: " << stp.read() << " sec" << endl;
}

//////////////////////////////////////////////////////////////////////////////

    } // namespace MutationUsage

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
