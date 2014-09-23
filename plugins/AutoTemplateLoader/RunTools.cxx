#include <iostream>

#ifdef CF_HAVE_BLITZ
  #include "RunTinyBlitz.hh"
  #include "RunBlitz.hh"
  #include "RunRoeFluxBlitz.hh"
  #include "RunRoeFluxTinyBlitz.hh"
#endif // CF_HAVE_BLITZ

#ifdef CF_HAVE_FLENS
  #include "RunFlens.hh"
  #include "RunRoeFluxFlens.hh"
#endif // CF_HAVE_FLENS

#ifdef CF_HAVE_TVMET
  #include "RunTVMET.hh"
  #include "RunRoeFluxTVMET.hh"
#endif // CF_HAVE_TVMET


#include "RunCFMathTools.hh"
#include "RunRoeFlux.hh"

using namespace std;

template < int ELEMS , int NTESTS >
double average_run(BaseRunLib * bench)
{
  double time = 0.;
  for (int i = 0; i < NTESTS; ++i)
  {
    time += bench->test(ELEMS);
  }
  return time / NTESTS;
}

template < int SIZE , int ELEMS , int NTESTS>
void runalgebra()
{
BaseRunLib * bench;

cout << SIZE << " " << ELEMS << " ";

#ifdef CF_HAVE_BLITZ
//   cout << "Blitz\t";
  bench = new RunBlitz<double,SIZE>();
  cout << average_run<ELEMS,NTESTS>(bench) << "\t";
  delete bench;
#endif

#ifdef CF_HAVE_FLENS
//   cout << "Flens\t";
  bench = new RunFlens<double,SIZE>();
  cout << average_run<ELEMS,NTESTS>(bench) << "\t";
  delete bench;
#endif

//   cout << "CF\t";
  bench = new RunCFMathTools<double,SIZE>();
  cout << average_run<ELEMS,NTESTS>(bench) << "\t";
  delete bench;

#ifdef CF_HAVE_BLITZ
//   cout << "TinyB\t";
  bench = new RunTinyBlitz<double,SIZE>();
  cout << average_run<ELEMS,NTESTS>(bench) << "\t";
  delete bench;
#endif

#ifdef CF_HAVE_TVMET
//   cout << "tvmet\t";
  bench = new RunTVMET<double,SIZE>();
  cout << average_run<ELEMS,NTESTS>(bench) << "\t";
  delete bench;
#endif

  cout << endl;
}

int main(int argc, char** argv)
{
BaseRunLib * bench;
// controls sensitivity to clock measure, more average, less sensitivity
#define AVERAGE 10
// controls sensitivity to lenght of tests, more meaningfull the bigger
#define NODES   800000

#if 0

  runalgebra<32,NODES,AVERAGE>();
  runalgebra<34,NODES,AVERAGE>();
  runalgebra<36,NODES,AVERAGE>();
  runalgebra<38,NODES,AVERAGE>();
  runalgebra<40,NODES,AVERAGE>();
  runalgebra<45,NODES,AVERAGE>();
  runalgebra<50,NODES,AVERAGE>();
#endif

#if 1

// must be 4
#define SIZE  4

  cout << NODES << "\t";

#ifdef CF_HAVE_BLITZ
//   cout << "Blitz\t";
  bench = new RunRoeFluxBlitz<double,SIZE>();
  cout << average_run<NODES,AVERAGE>(bench) << "\t";
  delete bench;
#endif

#ifdef CF_HAVE_FLENS
//   cout << "Flens\t";
  bench = new RunRoeFluxFlens<double,SIZE>();
  cout << average_run<NODES,AVERAGE>(bench) << "\t";
  delete bench;
#endif

//   cout << "CF\t";
  bench = new RunRoeFlux<double,SIZE>();
  cout << average_run<NODES,AVERAGE>(bench) << "\t";
  delete bench;

#ifdef CF_HAVE_BLITZ
//   cout << "TinyB\t";
  bench = new RunRoeFluxTinyBlitz<double,SIZE>();
  cout << average_run<NODES,AVERAGE>(bench) << "\t";
  delete bench;
#endif

#ifdef CF_HAVE_TVMET
//   cout << "tvmet\t";
  bench = new RunRoeFluxTVMET<double,SIZE>();
  cout << average_run<NODES,AVERAGE>(bench) << "\t";
  delete bench;
#endif

  cout << endl;

#endif

  return 0;
}

