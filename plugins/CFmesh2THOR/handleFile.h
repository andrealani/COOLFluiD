/////////////////////////////////////////////////////////////////////////////
/***************************************************************************/
/*                      handleFiles.h -  description                       */
/*                           -------------------                           */
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

#ifndef __HANDLEFILES_H__
#define __HANDLEFILES_H__

#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>

#include "smallCodeFunction.h"

FILE *opencfile(const char name[], char attribute);

void read_CFmesh_file();

void write_THOR_file();
void write_XYZ();

#endif
