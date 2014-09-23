
#include "MersenneTwister.hh"
#include <iostream>

//////////////////////////////////////////////////////////////////////////////

using namespace std;

//////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD {

  namespace RadiativeTransfer {

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

// initialization of static private members
unsigned long MersenneTwister::state[n] = {0x0UL};
int MersenneTwister::p = 0;
bool MersenneTwister::init = false;

//////////////////////////////////////////////////////////////////////////////

// generate new state vector
void MersenneTwister::New_State(){ 

  for (int i = 0; i < (n - m); ++i){
    state[i] = state[i + m] ^ recurrence(state[i], state[i + 1]);
  }

  for (int i = n - m; i < (n - 1); ++i){
    state[i] = state[i + m - n] ^ recurrence(state[i], state[i + 1]);
  }

  state[n - 1] = state[m - 1] ^ recurrence(state[n - 1], state[0]);
  p = 0; // reset position
}


//////////////////////////////////////////////////////////////////////////////

// init by array
void MersenneTwister::seed(const unsigned long* array, int size){ 

  seed(19650218UL);
  int i = 1, j = 0;
  for (int k = ((n > size) ? n : size); k; --k){
    state[i] = ((state[i] ^ ((state[i - 1] ^ (state[i - 1] >> 30)) * 1664525UL)) + array[j] + j) & 0xFFFFFFFFUL;
    ++j; j %= size;
    if ((++i) == n) { state[0] = state[n - 1]; i = 1; }
  }

  for (int k = n - 1; k; --k){
    state[i] = ((state[i] ^ ((state[i - 1] ^ (state[i - 1] >> 30)) * 1566083941UL)) - i) & 0xFFFFFFFFUL;
    if ((++i) == n) { state[0] = state[n - 1]; i = 1; }
  }

  state[0] = 0x80000000UL;
  p = n; // force New_State() to be called for next random number
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace RadiativeTransfer

} // namespace COOLFluiD




