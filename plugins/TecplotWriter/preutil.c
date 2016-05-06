#include "TECXXX.h"

#include <stdlib.h>    /* for malloc, mktemp */
#include <stdio.h>
#include <string.h>    /* strcmp, ... */
#include <sys/types.h>
#include <errno.h>
#include <unistd.h>    /* for unlink */

typedef enum   { GridMode,
                 WindowMode,
                 LogMode,
                 ScreenMode }                          CoordModeType;
typedef enum   { SolidLine,
                 DashedLine,
                 DashDotLine,
                 DottedLine,
                 LongDashLine,
                 DashDotDotLine }                      LineType;
typedef enum   { G_LineSegs,
                 G_Rectangle,
                 G_Circle,
                 G_Ellipse,
                 G_3DLineSegs}                         GeometryType;
typedef enum   { NoTextBox,
                 HollowTextBox,
                 FilledTextBox }                       TextBoxType;
typedef enum   { F_Helv,
                 F_HelvBold,
                 F_Greek,
                 F_Math,
                 F_UserDef,
                 F_Times,
                 F_TimesItalic,
                 F_TimesBold,
                 F_TimesItalicBold,
                 F_Courier,
                 F_CourierBold }                    FontType;
typedef enum   { S_Global,
                 S_Local }                             ScopeType;

typedef enum FormatType   {BLOCK,POINT,FEBLOCK,FEPOINT} FormatType;
typedef FILE     *FilePtr;

#define ZoneMarker        299.0
#define GeomMarker        399.0
#define TextMarker        499.0
#define CustomLabelMarker 599.0
#define EndHeaderMarker   357.0

#define MaxNumFiles    10
#define MaxChrsVarName 32
#define MaxChrsZneName 32

static LgIndex    DoDebug[MaxNumFiles] = {0,0,0,0,0,0,0,0,0,0};
static int        IsOpen[MaxNumFiles]  = {0,0,0,0,0,0,0,0,0,0};
static LgIndex    NumVars[MaxNumFiles];
static FNameType  DestFName[MaxNumFiles];
static FNameType  BlckFName[MaxNumFiles];
static FilePtr    BlckFile[MaxNumFiles];
static FilePtr    HeadFile[MaxNumFiles];
static FormatType ZoneFormat[MaxNumFiles];
static LgIndex    IMax[MaxNumFiles];
static LgIndex    JMax[MaxNumFiles];
static LgIndex    KMax[MaxNumFiles];
static int        NumIndicies[MaxNumFiles];
static LgIndex    DataCount[MaxNumFiles];
static int        CurZone[MaxNumFiles];
static int        CurFile = -1;


static void WriteI(F,I)
            FILE  *F;
            LgIndex *I;
{
  fwrite(I,4,1,F);
}


static void WriteR(F,R)
            FILE  *F;
            float *R;
{
  /* int    i; */
  /* char   bufr[4]; */

#ifndef CONVEX
  char  *buf = (char *)R;
#endif

#if defined CONVEX
  float  RR;
  char  *buf = (char *)&RR;

  RR = (float)dcvtid((double)*R);  /* Must use the Double precision
                                      version due to bug in single.
                                    */
#endif
  fwrite(buf,4,1,F);
}


static void WriteN(F,R)
            FILE  *F;
            float  R;
{
  float  RR;
  RR = R;
  WriteR(F,&RR);
}



static void WriteS(F,S)
            FILE  *F;
            char  *S;
{
  char *CPtr = S;
  while (*CPtr)
    {
      WriteN(F,(float)*CPtr);
      CPtr++;
    }
  WriteN(F,(float)0);
}





int  TECINI(
     char  *Title,
     char  *Variables,
     char  *FName,
     char  *ScratchDir,
     int   *Debug)
{
  int   I,L;
  char  RName[10];
  char *CPtr;
  int   IsOk = 1;
  int   NewFile = -1;

  for (I = 0; (I < MaxNumFiles) && (NewFile == -1); I++)
    {
      if (!IsOpen[I])
        NewFile = I;
    }

  if (NewFile == -1)
    {
      printf("Err: Too many files (%d) opened for printing\n",NewFile);
      return (-1);
    }

  if (CurFile == -1)
    CurFile = 0;

  DoDebug[NewFile] = *Debug;

  CurZone[NewFile] = 0;
  L = 0;
  if (FName != NULL)
    L = strlen(FName);
  if (L == 0)
    {
      printf("Err: Bad Filename for tecplot plot file, ref=%d\n",NewFile);
      return (-1);
    }
  DestFName[NewFile] = (char *)malloc(L+1);
  strcpy(DestFName[NewFile],FName);



  sprintf(RName,"tp%1dXXXXXX",NewFile+1);
  mktemp(RName);

  L = 0;
  if (ScratchDir != NULL)
    L = strlen(ScratchDir);
  BlckFName[NewFile] = (char *)malloc(L+1+8);
  if (ScratchDir != NULL)
    {
      strcpy(BlckFName[NewFile],ScratchDir);
      strcat(BlckFName[NewFile],"/");
    }
  else
    BlckFName[NewFile][0] = '\0';
  strcat(BlckFName[NewFile],RName);

  if (DoDebug[NewFile])
    {
      printf("Scratch File #%d: %s\n",NewFile+1,BlckFName[NewFile]);
      printf("Dest    File #%d: %s\n",NewFile+1,DestFName[NewFile]);
    }

  HeadFile[NewFile] = fopen(DestFName[NewFile],"w");
  BlckFile[NewFile] = fopen(BlckFName[NewFile],"w");

  if (BlckFile[NewFile] == NULL)
    {
      printf("Err: Cannot open scratch file for tecplot output....\n");
      printf("     Check permissions in scratch directory.\n");
      return (-1);
    }
  if (HeadFile[NewFile] == NULL)
    {
      printf("Err: Cannot open tecplot plot file... Check permissions.\n");
      return (-1);
    }


  WriteN(HeadFile[NewFile],6.3);

  WriteS(HeadFile[NewFile],Title);

  NumVars[NewFile] = 0;
  CPtr    = Variables;
  while (*CPtr)
    {
      while (*CPtr && ((*CPtr == ' ') || (*CPtr == ','))) CPtr++;
      if (*CPtr)
        {
          NumVars[NewFile]++;
          while (*CPtr && (*CPtr != ' ') && (*CPtr != ',')) CPtr++;
        }
    }

  if (NumVars[NewFile] == 0)
    {
      printf("Err: No variable names were defined\n");
      return (-1);
    }

  if (DoDebug[NewFile])
    printf("NumVars=%d\n",NumVars[NewFile]);

  WriteN(HeadFile[NewFile],(float)NumVars[NewFile]);
  CPtr    = Variables;
  while (*CPtr)
    {
      while (*CPtr && ((*CPtr == ' ') || (*CPtr == ','))) CPtr++;
      if (*CPtr)
        {
          while (*CPtr && (*CPtr != ' ') && (*CPtr != ','))
            {
              L = 0;
              if (L < MaxChrsVarName)
                WriteN(HeadFile[NewFile],(float)*CPtr);
              L++;
              CPtr++;
            }
          WriteN(HeadFile[NewFile],(float)0.0);
        }
    }
  IsOpen[NewFile] = 1;
  return (0);
}


static int CheckData()
{
  LgIndex NumDataPoints;

  if ((ZoneFormat[CurFile] == FEPOINT) ||
      (ZoneFormat[CurFile] == FEBLOCK))
    NumDataPoints = NumVars[CurFile]*IMax[CurFile];
  else
    NumDataPoints = NumVars[CurFile]*IMax[CurFile]*JMax[CurFile]*KMax[CurFile];

  if (NumDataPoints != DataCount[CurFile])
    {
      printf("Err in file %d: \n",CurFile+1);
      printf("     %d data values for Zone %d were processed \n",DataCount[CurFile],CurZone[CurFile]);
      printf("     %d data values were expected.\n",NumDataPoints);
      return (-1);
    }
  return (0);
}



int  TECZNE(
     char  *ZoneTitle,
     int   *IMx,
     int   *JMx,
     int   *KMx,
     char  *ZFormat)
{
  int        I;
  int        L = 0;
  int        IsOk = 1;
  char      *CPtr;

  if ((CurFile == -1) || (!IsOpen[CurFile]))
    {
      printf("Err:  TECZNE called for invalid file (%d)\n",CurFile+1);
      return (-1);
    }

  if (CurZone[CurFile] > 0)
    {
      if (CheckData() < 0)
        return (-1);
    }

  DataCount[CurFile] = 0;
  CurZone[CurFile]++;
  IMax[CurFile] = *IMx;
  JMax[CurFile] = *JMx;
  KMax[CurFile] = *KMx;

  WriteN(HeadFile[CurFile],(float)ZoneMarker);
  if (ZoneTitle != NULL)
    {
      WriteS(HeadFile[CurFile],ZoneTitle);
    }
  if (ZFormat == NULL)
    IsOk = 0;
  else if (!strcmp(ZFormat,"BLOCK"))
    ZoneFormat[CurFile] = BLOCK;
  else if (!strcmp(ZFormat,"FEBLOCK"))
    ZoneFormat[CurFile] = FEBLOCK;
  else if (!strcmp(ZFormat,"POINT"))
    ZoneFormat[CurFile] = POINT;
  else if (!strcmp(ZFormat,"FEPOINT"))
    ZoneFormat[CurFile] = FEPOINT;
  else
    IsOk = 0;
  if (IsOk == 0)
    {
      printf("Err: Bad Zone Format for Tecplot plot file %d\n",CurFile+1);
      return (-1);
    }

  if ((ZoneFormat[CurFile] == FEPOINT) ||
      (ZoneFormat[CurFile] == FEBLOCK))
    {
      if (KMax[CurFile] == 0)
        NumIndicies[CurFile] = 3;
      else if (KMax[CurFile] == 1)
        NumIndicies[CurFile] = 4;
      else if (KMax[CurFile] == 2)
        NumIndicies[CurFile] = 4;
      else if (KMax[CurFile] == 3)
        NumIndicies[CurFile] = 8;
    }


  WriteN(HeadFile[CurFile],(float)ZoneFormat[CurFile]);
  WriteN(HeadFile[CurFile],(float)-1.0);       /* No Zone Color Assignment */
  WriteN(HeadFile[CurFile],(float)IMax[CurFile]);
  WriteN(HeadFile[CurFile],(float)JMax[CurFile]);
  WriteN(HeadFile[CurFile],(float)KMax[CurFile]);

  WriteN(BlckFile[CurFile],(float)ZoneMarker);
  WriteN(BlckFile[CurFile],(float)0.0);

  if (DoDebug[CurFile])
    {
      printf("Writing Zone %d:\n",CurZone[CurFile]);
      printf("      Title = %s\n",ZoneTitle);
      printf("      Format= %s\n",ZFormat);
      printf("      IMax  = %d\n",IMax[CurFile]);
      printf("      JMax  = %d\n",JMax[CurFile]);
      printf("      KMax  = %d\n",KMax[CurFile]);
    }
  return (0);
}

static int CheckFile()
{
  if ((CurFile == -1) || (!IsOpen[CurFile]))
    {
      printf("Err:  Attempt to write tecplot data when no files are defined\n");
      return (-1);
    }
  return (1);
}


int TECDAT(
    LgIndex *N,
    float   *Data,
    long    *IsDouble)
{
  LgIndex I;

  if (CheckFile() == -1)
    return (-1);

  if (*IsDouble == 1)
    {
      double *d;
      float   f;
      d = (double *)Data;
      for (I = 0; I < *N; I++)
        {
          if (*d < -1.0E36)
            f = -1.0E36;
          else if ((*d > -1.0E-36) && (*d < 1.0E-36))
            f = 0.0;
          else
            f = *d;
          WriteR(BlckFile[CurFile],&f);
          d++;
        }
    }
  else
    {
      for (I = 0; I < *N; I++)
        WriteR(BlckFile[CurFile],&Data[I]);
    }
  DataCount[CurFile] += *N;
  return (0);
}

int TECNOD(
    LgIndex *NData)
{
  LgIndex L = NumIndicies[CurFile]*JMax[CurFile];
  LgIndex I;

  if ((CurFile == -1) || !IsOpen[CurFile])
    {
      printf("Err:  Attempt to write tecplot data on invalid file (%d)\n",CurFile+1);
      return (-1);
    }

  WriteN(BlckFile[CurFile],(float)0.0); /* No Repeat Variables */
  if ((ZoneFormat[CurFile] == BLOCK) || (ZoneFormat[CurFile] == POINT))
    {
      printf("Err: Tecplot output: Cannot call TECNOD if format is POINT or BLOCK\n");
      return (-1);
    }

  if (CheckData() < 0)
    return (-1);

  for (I = 0; I < L; I++)
    WriteI(BlckFile[CurFile],&NData[I]);
  return (0);
}



int  TECEND()
{
  char buf[4];

  if ((CurFile == -1) || (!IsOpen[CurFile]))
    {
      printf("Err:  Attempt to close invalid tecplot file %d\n",CurFile+1);
      return (-1);
    }

  if (CheckData() < 0)
    return (-1);

  WriteN(HeadFile[CurFile],(float)EndHeaderMarker);


  fclose(BlckFile[CurFile]);

  BlckFile[CurFile] = fopen(BlckFName[CurFile],"r");

/*
  rewind(BlckFile[CurFile]);
*/

  while (fread(buf,4,1,BlckFile[CurFile]))
    fwrite(buf,4,1,HeadFile[CurFile]);

  fclose(BlckFile[CurFile]);

  unlink(BlckFName[CurFile]);

  fclose(HeadFile[CurFile]);

  if (DoDebug[CurFile])
    printf("File %d closed.\n",CurFile+1);

  IsOpen[CurFile] = 0;
  CurFile = 0;
  while ((CurFile < MaxNumFiles) && !IsOpen[CurFile])
    CurFile++;

  if (CurFile == MaxNumFiles)
    CurFile = -1;

  return (0);
}




static void GetNextLabel(CPtr,NextLabel)
            char **CPtr;
            char *NextLabel;
{
  int N = 0;
  char *NPtr = NextLabel;
  *NPtr = '\0';
  /* Find label start */
  while ((**CPtr) && (**CPtr != '"'))
    (*CPtr)++;
  if (**CPtr)
    (*CPtr)++;
  while ((N < 60) && (**CPtr) && (**CPtr != '"'))
    {
      if (**CPtr == '\\')
        {
          (*CPtr)++;
        }
      *NPtr++ = **CPtr;
      N++;
      (*CPtr)++;
    }
  if (**CPtr)
    (*CPtr)++;
  *NPtr = '\0';
}


int TECLAB(char  *S)
{
  char     *CPtr = S;
  LgIndex   N = 0;
  char      Label[60];

  if (CheckFile() == -1)
    return (-1);

  if (DoDebug[CurFile])
    printf("\nInserting Custom Labels:\n");

  do
    {
      GetNextLabel(&CPtr,Label);
      if (*Label)
        N++;
    }
  while (*Label);

  if (N == 0)
    {
      printf("Err:  Invalid custom label string: %s\n",
             (S? S : " "));
      return (-1);
    }

  WriteN(HeadFile[CurFile],CustomLabelMarker);
  WriteN(HeadFile[CurFile],(float)N);
  CPtr = S;
  do
    {
      GetNextLabel(&CPtr,Label);
      if (*Label)
        {
          WriteS(HeadFile[CurFile],Label);
          if (DoDebug[CurFile])
            printf("          %s\n",Label);
        }
    }
  while (*Label);
  return (0);
}



static void WriteGeomDataBlock(Data,NumPts)
            float  *Data;
            LgIndex NumPts;
{
  LgIndex I;
  for (I = 0; I < NumPts; I++)
    WriteR(HeadFile[CurFile],&Data[I]);
}


static void ShowDebugColor(Color)
            LgIndex        Color;
{
  switch (Color)
    {
      case 0 : printf("BLACK\n"); break;
      case 1 : printf("RED\n"); break;
      case 2 : printf("GREEN\n"); break;
      case 3 : printf("BLUE\n"); break;
      case 4 : printf("CYAN\n"); break;
      case 5 : printf("YELLOW\n"); break;
      case 6 : printf("PURPLE\n"); break;
      case 7 : printf("WHITE\n"); break;
      case 15 :
      case 16 :
      case 17 :
      case 18 :
      case 19 :
      case 20 :
      case 21 :
      case 22 : printf("CUSTOM%1d\n",Color-14); break;
      default : printf("INVALID\n");
    }
}



int TECGEO(
    LgIndex *CoordMode,
    float   *XPos,
    float   *YPos,
    float   *ZPos,
    LgIndex *Scope,
    LgIndex *Zone,
    LgIndex *Color,
    LgIndex *FillColor,
    LgIndex *IsFilled,
    LgIndex *LineStyle,
    LgIndex *GeomType,
    LgIndex *NumSegments,
    LgIndex *NumSegPts,
    float   *XGeomData,
    float   *YGeomData,
    float   *ZGeomData)
{
  int I,S;

  if (CheckFile() == -1)
    return (-1);


  WriteN(HeadFile[CurFile],GeomMarker);
  WriteN(HeadFile[CurFile],(float)*CoordMode);
  WriteN(HeadFile[CurFile],(float)*Scope);
  WriteR(HeadFile[CurFile],XPos);
  WriteR(HeadFile[CurFile],YPos);
  WriteR(HeadFile[CurFile],ZPos);
  WriteN(HeadFile[CurFile],(float)*Zone);
  WriteN(HeadFile[CurFile],(float)*Color);
  WriteN(HeadFile[CurFile],(float)*FillColor);
  WriteN(HeadFile[CurFile],(float)*IsFilled);
  WriteN(HeadFile[CurFile],(float)*GeomType);
  WriteN(HeadFile[CurFile],(float)*LineStyle);
  if (DoDebug[CurFile])
    {
      printf("\nInserting Geometry\n");
      printf("  CoordMode    = ");
      switch ((CoordModeType)*CoordMode)
        {
          case GridMode   : printf("GRID\n"); break;
          case WindowMode : printf("WINDOW\n"); break;
          default : printf("INVALID\n");
        }

      printf("  Scope        = ");
      switch ((ScopeType)*Scope)
        {
          case S_Local    : printf("LOCAL\n"); break;
          case S_Global   : printf("GLOBAL\n"); break;
          default : printf("INVALID\n");
        }
      printf("  X,Y,Z Origin = %G,%G,%G\n",*XPos,*YPos,*ZPos);
      printf("  Zone         = %d\n",*Zone);
      printf("  Color        = "); ShowDebugColor(*Color);
      printf("  FillColor    = "); ShowDebugColor(*FillColor);
      printf("  IsFilled     = %d\n",*IsFilled);
      printf("  GeomType     = ");
      switch ((GeometryType)*GeomType)
        {
          case G_LineSegs    : printf("LINESEGS\n"); break;
          case G_Rectangle   : printf("RECTANGLE\n"); break;
          case G_Circle      : printf("CIRCLE\n"); break;
          case G_Ellipse     : printf("ELLIPSE\n"); break;
          case G_3DLineSegs  : printf("3D-LINESEGS\n"); break;
          default : printf("INVALID\n");
        }
      printf("  LineStyle    = ");
      switch ((LineType)*LineStyle)
        {
          case SolidLine      : printf("SOLID\n"); break;
          case DashedLine     : printf("DASHED\n"); break;
          case DashDotLine    : printf("DASHDOT\n"); break;
          case DottedLine     : printf("DOTTED\n"); break;
          case LongDashLine   : printf("LONGDASH\n"); break;
          case DashDotDotLine : printf("DASHDOTDOT\n"); break;
          default : printf("INVALID\n");
        }
    }
  if ((*GeomType == (int)G_LineSegs) ||
      (*GeomType == (int)G_3DLineSegs))
    {
      WriteN(HeadFile[CurFile],(float)*NumSegments);
      I = 0;
      for (S = 0; S < *NumSegments; S++)
        {
          WriteN(HeadFile[CurFile],(float)NumSegPts[S]);
          WriteGeomDataBlock(&XGeomData[I],NumSegPts[S]);
          WriteGeomDataBlock(&YGeomData[I],NumSegPts[S]);
          if (*GeomType == (int)G_3DLineSegs)
            WriteGeomDataBlock(&ZGeomData[I],NumSegPts[S]);
          I += NumSegPts[S];
        }
    }
  else if (*GeomType == (int)G_Circle)
    {
      WriteN(HeadFile[CurFile],XGeomData[0]);
    }
  else /* ellipse or rectangle */
    {
      WriteN(HeadFile[CurFile],XGeomData[0]);
      WriteN(HeadFile[CurFile],YGeomData[0]);
    }
  return (0);
}


int TECTXT(
    LgIndex *CoordMode,
    float   *XPos,
    float   *YPos,
    LgIndex *Scope,
    LgIndex *Font,
    float   *FontHeight,
    LgIndex *BoxType,
    float   *BoxMargin,
    LgIndex *BoxColor,
    LgIndex *BoxFillColor,
    LgIndex *Angle,
    LgIndex *Zone,
    LgIndex *TextColor,
    char    *Text)
{
  char *CPtr = Text;
  LgIndex L = strlen(Text);

  if (CheckFile() == -1)
    return (-1);

  WriteN(HeadFile[CurFile],TextMarker);
  WriteN(HeadFile[CurFile],(float)*CoordMode);
  WriteN(HeadFile[CurFile],(float)*Scope);
  WriteR(HeadFile[CurFile],XPos);
  WriteR(HeadFile[CurFile],YPos);
  WriteN(HeadFile[CurFile],(float)*Font);
  WriteR(HeadFile[CurFile],FontHeight);
  WriteN(HeadFile[CurFile],(float)*BoxType);
  WriteR(HeadFile[CurFile],BoxMargin);
  WriteN(HeadFile[CurFile],(float)*BoxColor);
  WriteN(HeadFile[CurFile],(float)*BoxFillColor);
  WriteN(HeadFile[CurFile],(float)*Angle);
  WriteN(HeadFile[CurFile],(float)*Zone);
  WriteN(HeadFile[CurFile],(float)*TextColor);
  WriteN(HeadFile[CurFile],(float)L);
  while (*CPtr)
    {
      WriteN(HeadFile[CurFile],(float)(*CPtr));
      CPtr++;
    }
  if (DoDebug[CurFile])
    {
      printf("\nInserting Text\n");
      printf("  CoordMode    = ");
      switch ((CoordModeType)*CoordMode)
        {
          case GridMode   : printf("GRID\n"); break;
          case WindowMode : printf("WINDOW\n"); break;
          default : printf("INVALID\n");
        }
      printf("  Scope        = ");
      switch ((ScopeType)*Scope)
        {
          case S_Local    : printf("LOCAL\n"); break;
          case S_Global   : printf("GLOBAL\n"); break;
          default : printf("INVALID\n");
        }
      printf("  X,Y   Origin = %G,%G\n",*XPos,*YPos);
      printf("  Font         = ");
      switch ((FontType)*Font)
        {
          case F_Helv            : printf("HELVETICA\n"); break;
          case F_HelvBold        : printf("HELVETICA-BOLD\n"); break;
          case F_Greek           : printf("GREEK\n"); break;
          case F_Math            : printf("MATH\n"); break;
          case F_UserDef         : printf("USERDEF\n"); break;
          case F_Times           : printf("TIMES\n"); break;
          case F_TimesItalic     : printf("TIMES-ITALIC\n"); break;
          case F_TimesBold       : printf("TIMES-BOLD\n"); break;
          case F_TimesItalicBold : printf("TIMES-ITALIC-BOLD\n"); break;
          case F_Courier         : printf("COURIER\n"); break;
          case F_CourierBold     : printf("COURIER-BOLD\n"); break;
          default : printf("INVALID\n");
        }
      printf("  Font Height  = %G\n",*FontHeight);
      printf("  Box Type     = ");
      switch ((TextBoxType)*BoxType)
        {
          case NoTextBox     : printf("NONE\n"); break;
          case HollowTextBox : printf("PLAIN\n"); break;
          case FilledTextBox : printf("FILLED\n"); break;
          default : printf("INVALID\n");
        }
      printf("  Box Margin   = %G\n",*BoxMargin);
      printf("  BoxColor     = "); ShowDebugColor(*BoxColor);
      printf("  BoxFillColor = "); ShowDebugColor(*BoxFillColor);
      printf("  Angle        = %d\n",*Angle);
      printf("  Zone         = %d\n",*Zone);
      printf("  TextColor    = "); ShowDebugColor(*TextColor);
      printf("  Text         = %s\n",Text);
    }
  return (0);
}


int  TECFIL( LgIndex   *F)
{
  if ((*F < 1) || (*F > MaxNumFiles))
    {
      printf("Err:  Invalid file number requested (%d)\n",(int) *F);
      return (-1);
    }
  if (!IsOpen[*F-1])
    {
      int I;
      printf("Err:  File set %d is not open!  File set not changed\n",(int) *F);
      printf("\n\nFile states are:\n");
      for (I = 0; I < MaxNumFiles; I++)
        printf("File %d, IsOpen=%d\n",I+1,IsOpen[I]);
      printf("Current File is: %d\n",CurFile+1);
      return (-1);
    }
  CurFile = *F-1;
  if (DoDebug[CurFile])
    {
      printf("Switching to tecplot file #%d\n\n",CurFile+1);
      printf("Current State is:\n");
      printf("  Debug     = %d\n",DoDebug[CurFile]);
      printf("  NumVars   = %d\n",NumVars[CurFile]);
      printf("  DestFName = %s\n",DestFName[CurFile]);
      printf("  BlckFName = %s\n",BlckFName[CurFile]);
      printf("  ZoneFormt = ");
      switch (ZoneFormat[CurFile])
        {
          case BLOCK    :
            {
              printf("BLOCK\n");
            } break;
          case POINT    :
            {
              printf("POINT\n");
            } break;
          case FEBLOCK  :
            {
              printf("FEBLOCK\n");
            } break;
          case FEPOINT  :
            {
              printf("FEPOINT\n");
            } break;
          default       :
            {
              printf("UNKNOWN\n");
            } break;
        }

      if ((ZoneFormat[CurFile] == BLOCK) ||
          (ZoneFormat[CurFile] == POINT))
        {
          printf("  IMax      = %d\n",(int) IMax[CurFile]);
          printf("  JMax      = %d\n",(int) JMax[CurFile]);
          printf("  KMax      = %d\n",(int) KMax[CurFile]);
        }
      else if ((ZoneFormat[CurFile] == FEBLOCK) ||
               (ZoneFormat[CurFile] == FEPOINT))
        {
          printf("  NumPoints = %d\n",(int) IMax[CurFile]);
          printf("  NumElmnts = %d\n",(int) JMax[CurFile]);
          printf("  EType     = ");
          if (KMax[CurFile] == 0)
            printf("Triangle\n");
          else if (KMax[CurFile] == 1)
            printf("Quadrilateral\n");
          else if (KMax[CurFile] == 2)
            printf("Tetrahedron\n");
          else if (KMax[CurFile] == 3)
            printf("Brick\n");
          else
            printf("Unknown\n");
        }
      printf("  DataCount = %d\n",DataCount[CurFile]);
      printf("  CurZone   = %d\n",CurZone[CurFile]);
    }
  return (0);
}

