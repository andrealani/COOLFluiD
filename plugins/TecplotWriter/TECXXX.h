#define DEC

#define TECINI  tecini
#define TECZNE  teczne
#define TECDAT  tecdat
#define TECNOD  tecnod
#define TECGEO  tecgeo
#define TECTXT  tectxt
#define TECLAB  teclab
#define TECFIL  tecfil
#define TECEND  tecend

#define round2(X)  ((X) >= 0 ? ((LgIndex)(X+0.49)) : ((LgIndex)(X-0.49)))

typedef long int LgIndex;

typedef char     *FNameType;


int  TECINI(
     char  *Title,
     char  *Variables,
     char  *FName,
     char  *ScratchDir,
     int   *Debug);

int  TECZNE(
     char  *ZoneTitle,
     int   *IMx,
     int   *JMx,
     int   *KMx,
     char  *ZFormat);

int TECDAT(
    LgIndex *N,
    float   *Data,
    long    *IsDouble);

int TECNOD(LgIndex *NData);

int  TECEND();

int TECLAB(char  *S);

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
    float   *ZGeomData);

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
    char    *Text);

int  TECFIL( LgIndex   *F);





