/////////////////////////////////////////////////////////////////////////////
/***************************************************************************/
/*                      CFmesh2Thor -  description                         */
/*                                                                         */
/*                       grid format convertor                             */
/*                  inspired by code of N. Villedieu                       */
/*                             -------------------                         */
/*  begin                : Tue Apr 29 15:11:29 CEST 2008                   */
/*  copyright            : (C) 2008 by Ing. Milan Zaloudek                 */
/*  email                : mili3@centrum.cz                                */
/***************************************************************************/
/////////////////////////////////////////////////////////////////////////////
/***************************************************************************/
/*                                                                         */
/*                      ALL RIGHTS RESERVED, (C) 2008                      */
/*                                                                         */
/***************************************************************************/
/////////////////////////////////////////////////////////////////////////////

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>

#include "handleFiles.h"
#include "smallCodeFunction.h"

/* GLOBAL VARIABLES */
uint DIM;
uint nbEqs;
uint nbNodes;
uint nbElem;
uint nbElemTypes;
uint nbElemPerType;
uint nbStates;
uint nbNodesPerType;
uint nbStatesPerType;
uint GeomPolyOrder;
uint SolPolyOrder;
bool SolutionPresent;

uint nbTRS;

ONE_NODE           *Node;
ONE_TETRA_ELEMENT  *Tetra;
ONE_TRS            *TRS;

void assign_boundary_to_element();

char inp_file[80];
char out_file[80];

int main(int argc, char *argv[])
{
if( argc == 3 ){
  strcpy(inp_file, argv[1]) ;
  strcpy(out_file, argv[2]) ;
//   strcpy(inp_file, "LightTipTailDrop5.CFmesh");
//   strcpy(out_file, "LightTipTailDrop5.thor0");

  cout << endl << "Input  file :  " << inp_file << endl;
  cout <<         "Output file :  " << out_file << endl << endl;

  cout << " --------------------------------" << endl
       << " ---------- reading -------------" << endl
       << " --------------------------------" << endl;
  read_CFmesh_file();



  cout << " --------------------------------" << endl
       << " ---- transforming to THOR ------" << endl
       << " --------------------------------" << endl;
  assign_boundary_to_element();



  cout << " --------------------------------" << endl
       << " ---------- writing -------------" << endl
       << " --------------------------------" << endl;
  write_THOR_file();
//   write_XYZ();
}

else{
  cout << "This program transforms CFmesh file (grid) to THOR file." << endl << endl
       << "  running command:   CFmesh2THOR <CFmesh file> <THOR file>" << endl
       << "                           CFmesh file must exist in actual folder" << endl
       << "                           THOR file will be created, containing transformed grid" << endl << endl
       << "  NOTE: several features (options) might not be implemented yet..." << endl
       << "               (C) written in Jul 2008 by Milan Zaloudek" << endl;
}
  return 0;

}
/******************************************************************************/

void assign_boundary_to_element(){

  cout << " Assigning boundary faces to elements ... " ;

  bool ElementFound;

  for( uint i=0; i<nbTRS; i++ ){
//     cout << endl << "  TRS - " << TRS[i].name ;
    for( uint j=0; j<TRS[i].size; j++ ){
      ElementFound = false;
      for( uint k=0; k<nbElem; k++ )
        if(ElementFound == true)  break;
        else{
          for( uint k1=0; k1<nbNodesPerType; k1++ )
            if( TRS[i].Face[j].NodeID[0] == Tetra[k].NodeID[k1] )

              for( uint k2=0; k2<nbNodesPerType; k2++ )
                if( TRS[i].Face[j].NodeID[1] == Tetra[k].NodeID[k2] )

                  for( uint k3=0; k3<nbNodesPerType; k3++ )
                    if( TRS[i].Face[j].NodeID[2] == Tetra[k].NodeID[k3] ){

                      TRS[i].Face[j].nb = Tetra[k].nb ;
                      ElementFound = true;
                    }
          }
    }
  }
  cout << "done." << endl << endl;
}
