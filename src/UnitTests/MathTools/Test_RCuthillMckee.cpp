#include <cstdio>
#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>

#include "Environment/CFEnv.hh"
#include "Common/ConnectivityTable.hh"
#include "Common/SwapEmpty.hh"

#include "MathTools/RCM.h"

using namespace std;
using namespace COOLFluiD;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Environment;







//******************************* MAIN ***************************************


int main ( int argc, char** argv )
{
  if ( argc != 2 )
  {
    std::cout << "Usage: " << argv[0] << " filename"  << std::endl;
    std::cout <<"\t Must provide a filename to open"  << std::endl;
    return 0;
  }

	CFEnv& cf_env = CFEnv::getInstance();  // build the environment
	cf_env.initiate ( argc, argv );        // initiate the environemnt
  Config::ConfigArgs args;
	cf_env.configure(args);
	cf_env.setup();


//............VARIABLES declaration........................
	ConnectivityTable<CFuint> cellnode;
	std::valarray < CFuint > new_id;
//.........................................................

	RCM::read_input ( argv[1], cellnode );     //CALL TO READ_INPUT

//.................checks .........................
// 	std::cout << cellnode;
//.................................................

	RCM::print_table ( "INPUT_numbering_TEC.dat", cellnode  ); //CALL TO PRINT_TABLE  ..to print the input

//.................checks .........................
// 	std::cout << cellnode;
// 	std::cout << "After printing input: " << std::endl;
//.................................................

	RCM::renumber ( cellnode, new_id );             //CALL TO RCM

  const CFuint nbnodes = new_id.size();
  for ( CFuint i=0; i< nbnodes; ++i ) {  std::cout<< new_id[i]<<endl;  }


        RCM::print_table ( "OUTPUT_numbering_TEC.dat", cellnode ); //CALL TO PRINT_TABLE  ..to print the output

//.................checks .........................
// 	std::cout << cellnode;
// 	std::cout << "After printing output: " << std::endl;
//.................................................

	return 0;
}



























