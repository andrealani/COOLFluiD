#! bashrc

clear

#echo "SCRIPT TO CONVERT THE ASCII INTO BINARIES"

#echo "Choose the files: (e.g.: /students/phd_ar/alaguna/LastIter/MRStudy2/T214/*.plt"
#read ASCIIFilesNames

#echo "Name of the output directory: (e.g.: ./T214)"
#read outputDirectory

#echo "Name of the output file: (e.g.: reconnection-time_214.plt) NOTE: Don't write the path"
#read outputName

#echo "Computing Postprocessing equations for reconnection case"

echo "SCRIPT TO CONVERT THE ASCII INTO BINARIES when the files are stored as T*/multi*.plt. 
      Limited for integer times"

for Directory in ./T*
do
  
  outputDirectory=$Directory
  ASCIIFilesNames=$Directory/multi*.plt
  time=$(echo ${Directory}| cut -c4-6)
  outputName="MRStudy4-time_"$time".plt"
  

  echo "Entering in $Directory"
  echo
  echo "Writing ./$Directory/$outputName"
  echo
  
  
  header="#!MC 1200"
  directoryVar="\$!VarSet |MFBD| = '"$outputDirectory"'"
  alterData="\$!ALTERDATA" 
  eq1="  EQUATION = '{Jz} = ({Ez} + {U0}*{By} - {V0}*{Bx})/0.075'"
  eq2="  EQUATION = '{J<sub>nonDim</sub>} = {Jz}/7.96e-3'"
  eq3="  EQUATION = '{n<sub>i</sub>}       = {rho0}/1.67262177774e-27'"
  eq4="  EQUATION = '{n<sub>n</sub>}       = {rho1}/1.6735327160314e-27'"
  eq5="  EQUATION = '{n<sub>Tot</sub>}    = {n<sub>i</sub>} + {n<sub>n</sub>} '"
  eq6="  EQUATION = '{rho<sub>Tot</sub>}  = {rho0} + {rho1}'"
  eq7="  EQUATION = '{Psi<sub>i</sub>}      = {n<sub>i</sub>} /({n<sub>Tot</sub>})'"
  eq8="  EQUATION = '{B}                           = sqrt({Bx}*{Bx} + {By}*{By} )'"
  eq9="  EQUATION = '{v<sub>A</sub>}       = {B}/(sqrt(4*3.14159265359*1e-7*{rho<sub>Tot</sub>}))'"
  eq10="  EQUATION = '{U0Ad} = {U0}/1.20071303e5'"
  eq11="  EQUATION = '{V0Ad} = {V0}/1.20071303e5'"
  eq12="  EQUATION = '{U1Ad} = {U1}/1.20071303e5'"
  eq13="  EQUATION = '{U1Ad} = {U1}/1.20071303e5'"
  eq14="  EQUATION = '{V1Ad} = {V1}/1.20071303e5'"
  outputVar="\$!WRITEDATASET \"|MFBD|/"$outputName"\""
  op1="INCLUDETEXT = NO"
  op2="INCLUDEGEOM = NO"
  op3="INCLUDECUSTOMLABELS = NO"
  op4="INCLUDEAUTOGENFACENEIGHBORS = YES"
  op5="ASSOCIATELAYOUTWITHDATAFILE = NO"
  op6="BINARY = YES"
  op7="USEPOINTFORMAT = NO"
  op8="PRECISION = 9"
  op9="TECPLOTVERSIONTOWRITE = TECPLOTCURRENT"
  end="\$!RemoveVar |MFBD|"


  echo $header > writeBinaryFile.mcr
  echo $directoryVar >> writeBinaryFile.mcr
  echo $alterData  >> writeBinaryFile.mcr
  echo $eq1 >> writeBinaryFile.mcr
  echo $alterData >> writeBinaryFile.mcr
  echo $eq2 >> writeBinaryFile.mcr
  echo $alterData >> writeBinaryFile.mcr
  echo $eq3 >> writeBinaryFile.mcr
  echo $alterData  >> writeBinaryFile.mcr
  echo $eq4 >> writeBinaryFile.mcr
  echo $alterData >> writeBinaryFile.mcr
  echo $eq5 >> writeBinaryFile.mcr
  echo $alterData >> writeBinaryFile.mcr
  echo $eq6 >> writeBinaryFile.mcr
  echo $alterData >> writeBinaryFile.mcr
  echo $eq7 >> writeBinaryFile.mcr
  echo $alterData >> writeBinaryFile.mcr
  echo $eq8 >> writeBinaryFile.mcr
  echo $alterData  >> writeBinaryFile.mcr
  echo $eq9 >> writeBinaryFile.mcr
  echo $alterData >> writeBinaryFile.mcr
  echo $eq10 >> writeBinaryFile.mcr
  echo $alterData >> writeBinaryFile.mcr
  echo $eq11 >> writeBinaryFile.mcr
  echo $alterData >> writeBinaryFile.mcr
  echo $eq12 >> writeBinaryFile.mcr
  echo $alterData  >> writeBinaryFile.mcr
  echo $eq13 >> writeBinaryFile.mcr
  echo $alterData >> writeBinaryFile.mcr
  echo $eq14 >> writeBinaryFile.mcr
  echo $outputVar >> writeBinaryFile.mcr
  echo $op1 >> writeBinaryFile.mcr
  echo $op2 >> writeBinaryFile.mcr
  echo $op3 >> writeBinaryFile.mcr
  echo $op4 >> writeBinaryFile.mcr
  echo $op5 >> writeBinaryFile.mcr
  echo $op6 >> writeBinaryFile.mcr
  echo $op7 >> writeBinaryFile.mcr
  echo $op8 >> writeBinaryFile.mcr
  echo $op9 >> writeBinaryFile.mcr
  echo $end >> writeBinaryFile.mcr

  echo "Finished writing the macro named ./writeBinaryFile.mcr"
  echo "Running tecplot"

  tec360 -b $ASCIIFilesNames -p writeBinaryFile.mcr 
done
echo "End of the script"
