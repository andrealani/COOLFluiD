#! bashrc

clear

echo "SCRIPT TO WRITE A PS from the plots"
echo
echo "NOTE: send the job from the directory where you have the ./T*"
echo
echo "Creating ./RESULTSplots where the plots are stored"
resultDirectory=./RESULTSplots
resultUion=$resultDirectory/Uion
resultVion=$resultDirectory/Vion
resultUneutral=$resultDirectory/Uneutral
resultVneutral=$resultDirectory/Vneutral
resultCurrent=$resultDirectory/Current
mkdir $resultDirectory
mkdir $resultUion $resultVion $resultUneutral $resultVneutral $resultCurrent

for Directory in ./T*
do
echo
echo "The name of the binary should be of the kind of MRStudy4-time_0.plt"
echo
  # Reading the time from the files and
  # storing it into variable $time
  inputDirectory= $Directory
  ASCIIFilesNames=$Directory/multi*.plt
  time=$(echo ${Directory}| cut -c4-6)
  inputFileName=$Directory/"MRStudy4-time_"$time".plt"
  

  echo "Entering in $Directory"
  echo
  echo "Writing the tecplot macro ./exportPS.mcr"
  echo
  echo "Generating plots for $inputFileName"

  ## Writing first part of the Macro
  header="#!MC 1200 "
  
  directoryVar="\$!VarSet |MFBD| = '"$inputDirectory"'"
  

  mirror="\$!CREATEMIRRORZONES 
  SOURCEZONES =  [1-256] 
  MIRRORVAR = 'Y' 
\$!ALTERDATA  [257-512] 
  EQUATION = '{V0} = (-1)*{V0}' 
\$!ALTERDATA  [257-512] 
  EQUATION = '{V1} = (-1)*{V1}' 
\$!ALTERDATA  [257-512] 
  EQUATION = '{Bx} = (-1)*{Bx}' 
\$!GLOBALCONTOUR 1  VAR = 3 
\$!CONTOURLEVELS RESETTONICE 
  CONTOURGROUP = 1 
  APPROXNUMVALUES = 15 
\$!FIELDLAYERS SHOWCONTOUR = YES 
\$!FIELDLAYERS SHOWSHADE = NO 
\$!FIELDLAYERS SHOWEDGE = NO 
\$!CREATEMIRRORZONES  
  SOURCEZONES =  [1-256] 
  MIRRORVAR = 'X' 
\$!CREATEMIRRORZONES  
  SOURCEZONES =  [257-512] 
  MIRRORVAR = 'X'
\$!ALTERDATA  [513-1024] 
  EQUATION = '{U0} = (-1)*{U0}' 
\$!ALTERDATA  [513-1024] 
  EQUATION = '{U1} = (-1)*{U1}' 
\$!ALTERDATA  [513-1024] 
  EQUATION = '{By} = (-1)*{By}' " 
  

  changeAxisAndCenter="\$!TWODAXIS DEPXTOYRATIO = 0.0250000000000000014 
\$!VIEW CENTER
\$!VIEW SETMAGNIFICATION 
  MAGNIFICATION = 6.17979502402 
\$!VIEW PUSH 
\$!TWODAXIS XDETAIL{TITLE{TITLEMODE = USETEXT}} 
\$!TWODAXIS XDETAIL{TITLE{TEXT = 'x [m]'}} 
\$!TWODAXIS YDETAIL{TITLE{TITLEMODE = USETEXT}} 
\$!TWODAXIS YDETAIL{TITLE{TEXT = 'y [m]'}} 
\$!RENAMEDATASETVAR  
  VAR = 20 
  NAME = 'Jz [A/m<sup>2</sup>]' 
\$!GLOBALCONTOUR 1  VAR = 20 
\$!CONTOURLEVELS RESETTONICE 
  CONTOURGROUP = 1 
  APPROXNUMVALUES = 15 "

  redraw="\$!REDRAW "
  

  setLegend="\$!GLOBALCONTOUR 1  LEGEND{SHOW = YES} 
\$!GLOBALCONTOUR 1  LEGEND{ISVERTICAL = NO} 
\$!GLOBALCONTOUR 1  LEGEND{AUTORESIZE = YES} 
\$!GLOBALCONTOUR 1  LEGEND{XYPOS{Y = 98}} "
  
  changeNameVars="\$!RENAMEDATASETVAR  
  VAR = 13 
  NAME = 'U<sub>ion</sub> [m/s]'
\$!RENAMEDATASETVAR  
  VAR = 14 
  NAME = 'V<sub>ion</sub> [m/s]' 
\$!RENAMEDATASETVAR  
  VAR = 15 
  NAME = 'U<sub>neutral</sub> [m/s]' 
\$!RENAMEDATASETVAR  
  VAR = 16 
  NAME = 'V<sub>neutral</sub> [m/s]' "


  setTime="\$!ATTACHTEXT  
  ANCHORPOS 
    { 
    X = 86.91489361702128 
    Y = 6.382978723404253 
    } 
  TEXT = 't = $time' 
\$!PICK SETMOUSEMODE 
  MOUSEMODE = SELECT 
\$!PICK ADDATPOSITION 
  X = 8.92718746066 
  Y = 7.67389525368 
  CONSIDERSTYLE = YES 
\$!PICK SHIFT 
  X = 0.0524235175626 
  Y = 0.209492635025 
\$!PICK ADDATPOSITION 
  X = 10.3688341936 
  Y = 7.73936170213 
  CONSIDERSTYLE = YES 
\$!PICK SETMOUSEMODE 
  MOUSEMODE = SELECT "

  #We start with the U_ion VAR = 13

  setContour="\$!GLOBALCONTOUR 1  VAR = 13 
\$!CONTOURLEVELS RESETTONICE 
  CONTOURGROUP = 1 
  APPROXNUMVALUES = 15 " 

  exportSetup="\$!EXPORTSETUP EXPORTFORMAT = PS 
\$!PRINTSETUP PALETTE = COLOR 
\$!EXPORTSETUP IMAGEWIDTH = 1275 "

  PSname="\$!EXPORTSETUP EXPORTFNAME = '$resultUion/Uions-time_$time.ps' "

  export="\$!EXPORT  
  EXPORTREGION = CURRENTFRAME"

  ending="\$!RemoveVar |MFBD| " 
  

  #Writing the macro
  echo "$header" > exportPS.mcr
  echo "$directoryVar" >> exportPS.mcr
  echo "$mirror" >> exportPS.mcr
  echo "$changeNameVars" >> exportPS.mcr
  echo "$changeAxisAndCenter" >> exportPS.mcr
  echo "$setTime" >> exportPS.mcr
  echo "$setContour" >> exportPS.mcr
  echo "$setLegend" >> exportPS.mcr
  echo "$redraw" >> exportPS.mcr
  echo "$exportSetup" >> exportPS.mcr
  echo "$PSname" >> exportPS.mcr
  echo "$export" >> exportPS.mcr

  #We continue the V_ion VAR = 14
  setContour="\$!GLOBALCONTOUR 1  VAR = 14 
\$!CONTOURLEVELS RESETTONICE 
  CONTOURGROUP = 1 
  APPROXNUMVALUES = 15 "
  PSname="\$!EXPORTSETUP EXPORTFNAME = '$resultVion/Vions-time_$time.ps' "

  echo "$setContour" >> exportPS.mcr  
  echo "$redraw" >> exportPS.mcr
  echo "$exportSetup" >> exportPS.mcr
  echo "$PSname" >> exportPS.mcr
  echo "$export" >> exportPS.mcr

  #We continue the U_neutral VAR = 15
  setContour="\$!GLOBALCONTOUR 1  VAR = 15 
\$!CONTOURLEVELS RESETTONICE 
  CONTOURGROUP = 1 
  APPROXNUMVALUES = 15 "
  PSname="\$!EXPORTSETUP EXPORTFNAME = '$resultUneutral/Uneutral-time_$time.ps' "

  echo "$setContour" >> exportPS.mcr  
  echo "$redraw" >> exportPS.mcr
  echo "$exportSetup" >> exportPS.mcr
  echo "$PSname" >> exportPS.mcr
  echo "$export" >> exportPS.mcr

  #We continue the V_neutral VAR = 16
  setContour="\$!GLOBALCONTOUR 1  VAR = 16 
\$!CONTOURLEVELS RESETTONICE 
  CONTOURGROUP = 1 
  APPROXNUMVALUES = 15 "
  PSname="\$!EXPORTSETUP EXPORTFNAME = '$resultVneutral/Vneutral-time_$time.ps' "

  echo "$setContour" >> exportPS.mcr  
  echo "$redraw" >> exportPS.mcr
  echo "$exportSetup" >> exportPS.mcr
  echo "$PSname" >> exportPS.mcr
  echo "$export" >> exportPS.mcr

  #We continue the Current VAR = 20
  setContour="\$!GLOBALCONTOUR 1  VAR = 20 
\$!CONTOURLEVELS RESETTONICE 
  CONTOURGROUP = 1 
  APPROXNUMVALUES = 15 "
  PSname="\$!EXPORTSETUP EXPORTFNAME = '$resultCurrent/Current-time_$time.ps' "

  echo "$setContour" >> exportPS.mcr  
  echo "$redraw" >> exportPS.mcr
  echo "$exportSetup" >> exportPS.mcr
  echo "$PSname" >> exportPS.mcr
  echo "$export" >> exportPS.mcr

  #ending the script
  echo "$ending" >> exportPS.mcr

  

  echo "Finished writing the macro named ./exportPS.mcr"
  echo "Running tecplot"

  tec360 -b $inputFileName -p exportPS.mcr 
done
echo "End of the script"
