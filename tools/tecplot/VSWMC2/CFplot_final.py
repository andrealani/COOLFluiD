def updateTecplotMacro():
   """
   updateTecplotMacro() : read the list of files and sort them by time
   """
   ofile = open("CFplot_final.mcr", "w")
   ifile = open("CFplot_template.mcr", "r")
   for line in ifile:
     if ("RESULTS_DIR" in line):
       newline = line.replace("RESULTS_DIR", sys.argv[1])			       
       ofile.write(newline)
     else:
       ofile.write(line)
     
   ifile.close()
   ofile.close()

#----------------------------------------------------------------------------------

def readFileList() -> list: 
   """
   readFileList() : read the list of files and sort them by time
   """
   fileList = []
   filepath = sys.argv[1] + "/list.txt"
   myfile = open(filepath, "r")
   beforeStr = sys.argv[2] + "-time_";
   afterStr  = ".plt";
   
   for line in myfile:
      ctime = line.replace(beforeStr,'')
      if (".surf" in line): ctime = ctime.replace(".surf",'')
      ctime = ctime.replace(afterStr,'')
      fname = line.replace('\n','')
      fileList.append((fname, int(ctime)))
      
   fileList.sort(key=lambda pair: pair[1])   
   return fileList

#----------------------------------------------------------------------------------

def runTecplotScript(fileList : list, surfList : list, cmd_changedir : str ):
   """
   runTecplotScript() : run the TECPLOT script to generate 
                        Bz_time_*.png  (Bz) 
                        p_time_*.png   (pressure) 
                        rho_time_*.png (density)
   """

   # update the TECPLOT script with the user-defined path
   updateTecplotMacro()
   
   cmd_tecplot = "module load tecplot/360_2009"
   system(cmd_tecplot)
   
   state = ["Bz", "rho", "p"] 
   counter = 0   
   for file in fileList:
      # generate output.plt.bin needed by TECPLOT script	        
      cmd_convert1 = cmd_tecplot + " ; " + cmd_changedir + " ; ln -sf " + file[0] + " output.plt ; preplot output.plt output.plt.bin"
      #print(cmd_convert1)
      system(cmd_convert1)

      # generate output.surf.plt.bin needed by TECPLOT script
      cmd_convert2 = cmd_tecplot + " ; " + cmd_changedir + " ; ln -sf " + surfList[counter][0] + " output.surf.plt ; preplot output.surf.plt output.surf.plt.bin"
      #print(cmd_convert2)
      counter+=1
      system(cmd_convert2)

      # generate the pics
      cmd_macro = cmd_tecplot + "; tec360 -b -p CFplot_final.mcr"
      print(cmd_macro)
      system(cmd_macro)

      for fvar in state:
       cmd_pics = cmd_tecplot + " ; " + cmd_changedir + " ; mv " + fvar + "_final.png" + " " + fvar + "_time_" + str(file[1]) + ".png"
       print(cmd_pics)
       system(cmd_pics)

#gimp *.png&

#----------------------------------------------------------------------------------

# python3 CFplot_final.py /Users/lani/Projects/VSWMC2/COOLFluiD/RESULTS Unsteady06042000Storm_Final
if __name__ == "__main__":
   from os import system
   import sys

   # dump the list of TECPLOT files 
   cmd_files = "ls " + sys.argv[2] + "-time_*.plt >& list.txt";
   cmd_changedir = "cd " + sys.argv[1]
   cmd_path = cmd_changedir + " ; $PWD ; " +  cmd_files;
   system(cmd_path)
   
   # read the list of files and sort them
   fileList = readFileList()

   # separate the surface files from the volume files
   surfList = []
   for file in fileList:
      if ("surf" in file[0]): surfList.append(file) ; fileList.remove(file) 
			    
   print(".surf.plt files: ", surfList)
   print(".plt      files: ", fileList)
   
   if (len(fileList) != len(surfList)):
      print("\nERROR: Number of *.plt [", len(fileList) , "] != number of *.surf.plt [" , len(surfList), "]\n")
      exit(1)

   runTecplotScript(fileList, surfList, cmd_changedir)
    
