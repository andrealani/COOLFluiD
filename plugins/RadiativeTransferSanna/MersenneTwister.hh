#ifndef COOLFluiD_RadiativeTransfer_MersenneTwister_hh
#define COOLFluiD_RadiativeTransfer_MersenneTwister_hh

//////////////////////////////////////////////////////////////////////////////

#include <iostream>

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

 namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

  /**
   * This class implement the Mersenne-Twister MT19937 random number generator
   * Based on code by Makoto Matsumoto and Takuji Nishimura.
   * The Mersenne Twister is an algorithm for generating random numbers,
   * designed with consideration of the flaws in various other generators.
   * It a very fast random number generator, with a maximum period of 2^19937-1
   * and an order of equidistribution of 623 dimensions. 
   * The generator is also fast; it avoids multiplication and
   * division, and it benefits from caches and pipelines.  
   * For more information see the inventors' web page at
   * http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html

   * Reference
   * M. Matsumoto and T. Nishimura, "Mersenne Twister: A 623-Dimensionally
   * Equidistributed Uniform Pseudo-Random Number Generator", ACM Transactions on
   * Modeling and Computer Simulation, Vol. 8, No. 1, January 1998, pp 3-30.
   *
   * @author Alessandro Sanna
   *
   */
class MersenneTwister                          
{
public:

  /**
   * Default constructor: uses default seed only if this is the first instance
   */
  MersenneTwister() : 
    m_c0(1.0 / 4294967295.0), 
    m_c1(1.0 / 4294967296.0), 
    m_c2(1.0 / 9007199254740992.0)
  { 
    if (!init) seed(5489UL); 
    init = true;
  }

  /**
   * Constructor: with 32 bit int as seed
   */
  MersenneTwister(unsigned long s) :
    m_c0(1.0 / 4294967295.0), 
    m_c1(1.0 / 4294967296.0), 
    m_c2(1.0 / 9007199254740992.0)
  {
    seed(s); 
    init = true;
  }
  
  /**
   * Constructor: with array of size 32 bit ints as seed
   */
  MersenneTwister(const unsigned long* array, int size) :
    m_c0(1.0 / 4294967295.0), 
    m_c1(1.0 / 4294967296.0), 
    m_c2(1.0 / 9007199254740992.0)
  {
    seed(array, size); 
    init = true;
  }
  
  /**
   * Default destructor
   */
   virtual ~MersenneTwister() {}

  /**
   * Seed with 32 bit integer
   */
   void seed(unsigned long);

  /**
   * Seed with array
   */
   void seed(const unsigned long*, int size); 

  /**
   * generate 32 bit random integer
   */
   unsigned long rand();

  /**
   * generates double floating point numbers in the closed interval [0, 1]
   */
   double rand_c0c1();

  /**
   * generates double floating point numbers in the closed interval [0, n]
   */
   double rand_c0cn( const double n );

  /**
   * generates double floating point numbers in the half-open interval [0, 1)
   */
   double rand_c0o1();
  
  /**
   * generates double floating point numbers in the half-open interval [0, n)
   */
   double rand_c0on( const double n );

  /**
   * generates double floating point numbers in the open interval (0, 1)
   */
   double rand_o0o1();

  /**
   * generates double floating point numbers in the open interval (0, n)
   */
   double rand_o0on( const double n );

  /**
   * generates 53 bit resolution doubles in the half-open interval [0, 1)
   */
   double rand_53();

  /**
   * Overload operator() to make this a generator (functor)
   */
   double operator()() { return rand_c0c1(); }

  
private: 
  
  // constant value (2^32 - 1)
  const double m_c0; 
  
  // constant value (2^32)
  const double m_c1;
  
  // constant value (from Isaku Wada)
  const double m_c2;
  
  // compile time constants
  static const int n = 624, m = 397;
  
  // the variables below are static (no duplicates can exist):
  
  // state vector array
  static unsigned long state[n];
  
    // position in state array
  static int p;
  
  // true if init function is called
  static bool init;
  
  // private functions used to generate the pseudo random numbers:
  
  // used by gen_state()
  unsigned long recurrence(unsigned long, unsigned long);
  
  // generate new state
  void New_State();
  
  // make copy constructor and assignment operator unavailable, they don't make sense:
  
  // copy constructor not defined
  MersenneTwister(const MersenneTwister&);
  
  // assignment operator not defined
  void operator=(const MersenneTwister&);
  
}; // end of class MersenneTwister

//////////////////////////////////////////////////////////////////////////////
// inline for speed, must therefore reside in header file

inline unsigned long MersenneTwister::recurrence(unsigned long u, unsigned long v) {
  return (((u & 0x80000000UL) | (v & 0x7FFFFFFFUL)) >> 1)
    ^ ((v & 1UL) ? 0x9908B0DFUL : 0x0UL);
}

inline unsigned long MersenneTwister::rand() { // generate 32 bit random int
  if (p == n /*new state vector needed*/){
    New_State();
    //std::cout<<"\n new state \n"<<std::endl;
  } 
// New_State() is split off to be non-inline, because it is only called once
// in every 624 calls and otherwise rand() would become too big to get inlined
// below: tempering in order to increase the equidistribution
  unsigned long x = state[p++];
  x ^= (x >> 11);
  x ^= (x << 7) & 0x9D2C5680UL;
  x ^= (x << 15) & 0xEFC60000UL;
  return x ^ (x >> 18);
}

// generates double floating point numbers in the closed interval [0, 1]
inline double MersenneTwister::rand_c0c1()
{ return static_cast<double>(rand()) *m_c0;} // divided by 2^32 - 1
   
// generates double floating point numbers in the closed interval [0, n]
inline double MersenneTwister::rand_c0cn( const double n )
{ return rand_c0c1() * n; }

// generates double floating point numbers in the half-open interval [0, 1)
inline double MersenneTwister::rand_c0o1()
{ return static_cast<double>(rand()) *m_c1; } // divided by 2^32
   
// generates double floating point numbers in the half-open interval [0, n)
inline double MersenneTwister::rand_c0on( const double n )
{ return rand_c0o1() * n; }

// generates double floating point numbers in the open interval (0, 1)
inline double MersenneTwister::rand_o0o1()
{ return (static_cast<double>(rand()) + 0.5) * m_c1; } // divided by 2^32
   
// generates double floating point numbers in the open interval (0, n)
inline double MersenneTwister::rand_o0on( const double n )
{ return rand_o0o1() * n; }

// generates 53 bit resolution doubles in the half-open interval [0, 1)
inline double MersenneTwister::rand_53()
{
  unsigned long a = static_cast<double>(rand() >> 5);
  unsigned long b = static_cast<double>(rand() >> 6);
  return ( a * 67108864.0 + b ) * m_c2;  // by Isaku Wada
}
   
// init by 32 bit seed
inline void MersenneTwister::seed(unsigned long s)
{  
  state[0] = s & 0xFFFFFFFFUL;
  for (int i = 1; i < n; ++i){
    state[i] = (1812433253UL * (state[i - 1] ^ (state[i - 1] >> 30)) + i) & 0xFFFFFFFFUL;
  }
  
  p = n; // force New_State() to be called for next random number
}
   
//////////////////////////////////////////////////////////////////////////////

  } // namespace RadiativeTransfer

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_RadiativeTransfer_MersenneTwister_hh
