
///@author Milan Zaloudek

#include "handleFiles.h"

FILE *opencfile(const char name[30], char attribute){
  FILE *ident;
  switch (attribute){
      case 'w' : ident = fopen(name , "w");break;
      case 'r' : ident = fopen(name , "r");break;
      case 'a' : ident = fopen(name , "a");break;
      default  : printf("\nInvalid attribute \"%c\" was found, while opening \" %s \"\n",attribute,name);exit(1);break;
  }
  if(!ident){printf("\n\nError occured during \" %s \" opening.\n",name);exit(1);}
  return ident;
}

void read_CFmesh_header(FILE *input){

  char  dummy[50];
  int   dummyInt;
  float dummyFloat;
  char  dummyString[50];

  extern uint DIM;
  extern uint nbEqs;
  extern uint nbNodes;
  extern uint nbStates;
  extern uint nbElem;
  extern uint nbElemTypes;
  extern uint nbElemPerType;
  extern uint GeomPolyOrder;
  extern uint SolPolyOrder;
  extern uint nbNodesPerType;
  extern uint nbStatesPerType;
  extern bool SolutionPresent;

  do{ fscanf(input, "%s", &dummy) ;

      if( strcmp(dummy,"!CFMESH_FORMAT_VERSION")==0 ){
        fscanf(input, "%f", &dummyFloat);
        cout << "  CFmesh format :  " << dummyFloat << endl;
        if( fabs(dummyFloat-1.3) > 0.001 ) cout << endl << "  This convertor is optimized for version 1.3" << endl
                                             << "     some functions couldn't work properly!" << endl << endl;
      }

      if( strcmp(dummy,"!COOLFLUID_VERSION")==0 ){
        fscanf(input, "%s", &dummyString);
        cout << "  obtained by CF ver. " << dummyString << endl;
      }

      if( strcmp(dummy,"!COOLFLUID_SVNVERSION")==0 ){
        fscanf(input, "%s", &dummyString);
        cout << "  obtained by CF ver. " << dummyString << "(svn)" << endl;
      }

      if( strcmp(dummy,"!NB_DIM")==0 ){
        fscanf(input, "%d", &DIM);
        cout << "  Dimensions :  " << DIM << endl;
      }

      if( strcmp(dummy,"!NB_EQ")==0 ){
        fscanf(input, "%d", &nbEqs);
        cout << "  Number of equations :  " << nbEqs << endl;
      }

      if( strcmp(dummy,"!NB_NODES")==0 ){
        fscanf(input, "%d", &nbNodes);
        fscanf(input, "%d", &dummyInt);
        cout << "  Number of nodes :  " << nbNodes << endl;
      }

      if( strcmp(dummy,"!NB_STATES")==0 ){
        fscanf(input, "%d", &nbStates);
        fscanf(input, "%d", &dummyInt);
        cout << "  Number of states :  " << nbStates << endl;
      }

      if( strcmp(dummy,"!NB_ELEM")==0 ){
        fscanf(input, "%d", &nbElem);
        cout << "  Number of elements :  " << nbElem << endl;
      }

      if( strcmp(dummy,"!NB_ELEM_TYPES")==0 ){
        fscanf(input, "%d", &nbElemTypes);
        cout << "  Number of element types present :  " << nbElemTypes << endl;
      }

      if( strcmp(dummy,"!GEOM_POLYORDER")==0 ){
        fscanf(input, "%d", &GeomPolyOrder);
        cout << "  Polynomial order of elements :  " << GeomPolyOrder << endl;
      }

      if( strcmp(dummy,"!SOL_POLYORDER")==0 ){
        fscanf(input, "%d", &SolPolyOrder);
        cout << "  Polynomial order of solution :  " << SolPolyOrder << endl;
      }

      if( strcmp(dummy,"!ELEM_TYPES")==0 ){
        fscanf(input, "%s", &dummy);
        if(strcmp(dummy,"Tetra")!=0){
          cout << "  Unhandled element type. Program stops. Besides - THOR can't handle anything but tetras" << endl;
        }
        else{
          cout << "  Elements :  Tetrahedrons" << endl;
        }
      }

      if( strcmp(dummy,"!NB_ELEM_PER_TYPE")==0 ){
        fscanf(input, "%d", &nbElemPerType);
        cout << "  Number of elements per type :  " << nbElemPerType << endl;
      }

      if( strcmp(dummy,"!NB_NODES_PER_TYPE")==0 ){
        fscanf(input, "%d", &nbNodesPerType);
        cout << "  Number of nodes per types :  " << nbNodesPerType << endl;
      }

      if( strcmp(dummy,"!NB_STATES_PER_TYPE")==0 ){
        fscanf(input, "%d", &nbStatesPerType);
        cout << "  Number of states per types :  " << nbStatesPerType << endl;
      }

      if( strcmp(dummy,"!LIST_STATE")==0 ){
        fscanf(input, "%d", &dummyInt);
        if( dummyInt == 0 ){ 
              SolutionPresent = false;
        }
        else{ SolutionPresent = true;
              cout << "  CFmesh file contains solution." << endl;
        }
      }

    }while( strcmp(dummy,"!LIST_ELEM")!=0 );
}

///////////////////////////////////////////////////////////////////////////////

void read_CFmesh_file()
{

  extern char inp_file[];

  cout << endl << " Reading CFmesh file." << endl;
  FILE *input = opencfile(inp_file,'r');


  char dummy[50];
  int dummyNumber;

//   extern uint DIM;
//   extern uint nbEqs;
  extern uint nbNodes;
  extern uint nbStates;
  extern uint nbElem;
  extern uint nbElemTypes;
  extern uint nbElemPerType;
  extern uint GeomPolyOrder;
  extern uint SolPolyOrder;
  extern uint nbNodesPerType;
  extern uint nbStatesPerType;
  extern bool SolutionPresent;

  /* Read header *****************************/
  cout << endl << "  Reading header" << endl;

      read_CFmesh_header(input);


  /* Read elements ***************************/
  cout << endl << " Reading elements ... " ;
  extern ONE_TETRA_ELEMENT *Tetra;
  Tetra = allocateTetra(nbElem);

  for( uint i=0; i<nbElem; i++){

    Tetra[i].nb = i+1;

    for( uint j=0; j<nbNodesPerType; j++){
      fscanf(input,"%d", &Tetra[i].NodeID[j]);
      Tetra[i].NodeID[j] += 1;  // THOR starts numbering with '1' while CFmesh with '0'
    }
    for( uint j=0; j<nbStatesPerType; j++){
      fscanf(input,"%d", &Tetra[i].StateID[j]);
      Tetra[i].StateID[j] += 1; // THOR starts numbering with '1' while CFmesh with '0'
    }
  }
  cout << "done." << endl << endl;

  /* Read boundaries *************************/
  extern uint nbTRS;
  extern ONE_TRS *TRS;

  fscanf(input, "%s", &dummy) ;

  if( strcmp(dummy,"!NB_TRSs")==0 ){
    fscanf(input, "%d", &nbTRS);
    cout << " Number of boundary segments :  " << nbTRS << endl;
    TRS = allocateTRS( nbTRS );
  }

  for( uint i=0; i<nbTRS; i++ ){
    do{
      fscanf(input, "%s", &dummy) ;

      if( strcmp(dummy,"!TRS_NAME")==0 ){
        fscanf(input, "%s", &TRS[i].name);
        cout << "  Loading part :  " << TRS[i].name << endl;
      }

      if( strcmp(dummy,"!NB_TRs")==0 )
        fscanf(input, "%i", &TRS[i].nb);

      if( strcmp(dummy,"!NB_GEOM_ENTS")==0 )
        fscanf(input, "%i", &TRS[i].size);

      if( strcmp(dummy,"!GEOM_TYPE")==0 )
        fscanf(input, "%s", &TRS[i].type);

    }while( strcmp(dummy,"!LIST_GEOM_ENT")!=0 );

    uint NodeSize, StateSize;

    for( uint j=0; j<TRS[i].size; j++ ){
      fscanf(input, "%d %d", &NodeSize, &StateSize);
      for( uint k=0; k<NodeSize; k++){
        fscanf(input, "%d", &TRS[i].Face[j].NodeID[k]);
        TRS[i].Face[j].NodeID[k] += 1; // THOR starts numbering with '1' while CFmesh with '0'
      }
      for( uint k=0; k<StateSize; k++){
        fscanf(input, "%d", &TRS[i].Face[j].StateID[k]);
        TRS[i].Face[j].StateID[k] += 1; // THOR starts numbering with '1' while CFmesh with '0'
      }
    }
  }

  fscanf(input, "%s", &dummy) ;
  if( strcmp(dummy,"!LIST_NODE")!=0 ){
    cout << endl << "  Wrong TRS loading - finished either to early or to late..." << endl ;
    exit(1);
  }
  /* Read node coordinates *******************/
  extern ONE_NODE *Node;

  Node = allocateCell(nbNodes);

  cout << " Reading node coordinates ... " ;
  for( uint i=0; i<nbNodes; i++){
    fscanf(input,"%lf %lf %lf", &Node[i].x, &Node[i].y, &Node[i].z);
  }
  cout << "done." << endl;

  //---------------------------------------------------------------------------

  fscanf(input, "%s", &dummy) ;

  if( strcmp(dummy,"!LIST_STATE")==0 ){
    fscanf(input, "%d", &dummyNumber);
    if( dummyNumber == 0 ){ 
          SolutionPresent = false;
          cout << "  No solution present." << endl;
    }
    else{ SolutionPresent = true;
          cout << "  CFmesh file contains solution." << endl;
    }
  }

  fscanf(input, "%s", &dummy) ;
  if( strcmp(dummy,"!END")!=0 ){
    cout << endl << " Loading not complete..." << endl << endl;
    exit(1);
  }
  else{
    cout << endl << " CFmesh file successfully loaded." << endl << endl;
  }

  fclose(input);

}

///////////////////////////////////////////////////////////////////////////////

void write_THOR_file(){

  extern char out_file[];

  cout << " Assembling THOR file." << endl;
  FILE *output = opencfile(out_file,'w');

  cout << "  writing header ... " ;
  extern uint DIM;
  extern uint nbElem;
  extern uint nbNodes;
  extern uint nbTRS;
  extern ONE_TRS *TRS;
  extern bool SolutionPresent;

  fprintf(output, "%d", DIM );
  if( SolutionPresent == false ){
    fprintf(output, " 0 0");
  }
  else{
    cout << " Unhandled exception - SolutionPresent\n Program stops." << endl;
    exit(0);
  }

  fprintf(output, "\n%d", nbElem );
  fprintf(output, " %d", nbNodes );

  uint sum = 0;
  for( uint i=0; i<nbTRS; i++ )
    sum += TRS[i].size;

  fprintf(output, " %d", sum );
  fprintf(output, " %d", nbTRS );
  cout << "done." << endl;

  extern uint nbElemTypes;
  extern uint nbNodesPerType;
  extern ONE_TETRA_ELEMENT *Tetra;

  cout << "  writing connectivity ... " ;
  fprintf(output, "\n%d", nbElemTypes );
  fprintf(output, "\n%d %d", nbNodesPerType, nbElem );
  for( uint i=0; i<nbElem; i++ ){
    fprintf(output, "\n");
    for( uint j=0; j<nbNodesPerType; j++){
      fprintf(output, "%d ", Tetra[i].NodeID[j] );
    }
  }

  for( uint i=0; i<nbTRS; i++ ){
    const uint TRScode = i+1;
    fprintf(output, "\n%d %d", TRScode, TRS[i].size );
    for( uint j=0; j<TRS[i].size; j++ ){
      fprintf(output, "\n%d 3 %d %d %d", TRS[i].Face[j].nb, TRS[i].Face[j].NodeID[0], TRS[i].Face[j].NodeID[1], TRS[i].Face[j].NodeID[2] );
    }
  }
  cout << "done." << endl;

  //-------------------------------------------------------

  cout << "  writing nodal data ... " ;

  extern ONE_NODE *Node;

  fprintf(output, "\n");
  for( uint i=1; i<=nbNodes; i++ ){
    fprintf(output, "%d ", i);
    if( (i%10)==0 ) fprintf(output, "\n");
  }

  fprintf(output, "\n");
  for( uint i=0; i<nbNodes; i++ ){
    fprintf(output, "%15.12lf ", Node[i].x);
    if( (i%4)==0 ) fprintf(output, "\n");
  }
  for( uint i=0; i<nbNodes; i++ ){
    fprintf(output, "%15.12lf ", Node[i].y);
    if( (i%4)==0 ) fprintf(output, "\n");
  }
  for( uint i=0; i<nbNodes; i++ ){
    fprintf(output, "%15.12lf ", Node[i].z);
    if( (i%4)==0 ) fprintf(output, "\n");
  }

  cout << "done." << endl << endl;

  fclose(output);
}

///////////////////////////////////////////////////////////////////////////////

void write_XYZ(){

  extern char out_file[];

  FILE *output = opencfile(out_file,'w');

  extern uint nbNodes;
  extern ONE_NODE *Node;

  for( uint i=0; i<nbNodes; i++ )
    fprintf(output, "\n%lf  %lf  %lf", Node[i].x, Node[i].y, Node[i].z);

  fclose(output);
}
