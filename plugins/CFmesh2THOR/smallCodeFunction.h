/////////////////////////////////////////////////////////////////////////////
/***************************************************************************/
/*                 smallCodeFunction.cpp  -  description                   */
/*                                                                         */
/*               external functions with programming basis                 */
/*                                                                         */
/*                        -------------------                              */
/*  begin                : Tue Jul 08 14:16:29 CEST 2008                   */
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

#ifndef __SCODEFUNCTION_H__
#define __SCODEFUNCTION_H__

#include <iostream.h>
#include <math.h>

#include "handleFiles.h"


typedef struct{
  double x;
  double y;
  double z;
}ONE_NODE;


typedef struct{
  uint nb;
  uint NodeID[4];
  uint StateID[4];
}ONE_TETRA_ELEMENT;


typedef struct{
  char name[20];
  char type[20];
  uint nb;
  uint size;
  ONE_TETRA_ELEMENT Face[15000];
}ONE_TRS;

double * allocateVector(int n);

double **allocateMatrix(int n, int m);

double ***allocateTensor(int n, int m, int k);

ONE_NODE *allocateCell(int nb);

ONE_TETRA_ELEMENT *allocateTetra(int nb);

ONE_TRS *allocateTRS(int nb);

#endif
