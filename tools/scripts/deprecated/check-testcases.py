#!/usr/bin/env python

import os                       # operating system operations
import string                   # string operations
import sys                      # system operations
import commands                 # for running programs
import time                     # for timing
import datetime                 # for date

class StopWatch:
    """ Implements a stop watch."""
    def __init__(self, parent=None):
        self._start = 0.0
        self._elapsedtime = 0.0
        self._running = 0
        self.timestr = ""

    def getStr(self):
        """ Print the time label. """
        self._setTime(self._elapsedtime)
        return self.timestr

    def _update(self):
        """ Update the label with elapsed time. """
        self._elapsedtime = time.time() - self._start
        self._setTime(self._elapsedtime)

    def _setTime(self, elap):
        """ Set the time string to Minutes:Seconds:Hundreths """
        minutes = int(elap/60)
        seconds = int(elap - minutes*60.0)
        hseconds = int((elap - minutes*60.0 - seconds)*100)
        self.timestr = '%02d:%02d:%02d' % (minutes, seconds, hseconds)

    def Start(self):
        """ Start the stopwatch, ignore if running. """
        if not self._running:
            self._start = time.time() - self._elapsedtime
            self._update()
            self._running = 1

    def Stop(self):
        """ Stop the stopwatch, ignore if stopped. """
        if self._running:
            self._elapsedtime = time.time() - self._start
            self._setTime(self._elapsedtime)
            self._running = 0

    def Reset(self):
        """ Reset the stopwatch. """
        self._start = time.time()
        self._elapsedtime = 0.0
        self._setTime(self._elapsedtime)

# global variables
fuzz = 10E-4                     # the accuracy for the residual comparison
convname = 'convergence.plt'     # the convergence file to open
PROG="./run-coolfluid.sh"        # the solver program/script to run
OUTFILE=".testcases.output.tmp"  # the temporary file to put the output of each run
CATALOGFILE="testcases.catalog"  # the file with the list of testcases
LOGFILE="testcases.log"          # the file to put the log of all the runs
ERRORFILE="testcases.err"        # the file to put the log of all the runs that FAIL
CHECKLOG="check-testcases.log"   # file to put the summary log
fout = open(CHECKLOG, 'w')

def print_stdout(msg):
	sys.stdout.write(msg)
	fout.write(msg)

def print_output(case,stpwatch,output,residual):
  """function to do the output"""
  smsg = ''
  time = stpwatch.getStr()
  lmsg = ''
  if(output[0] != 0):
    smsg = 'FAIL'
    lmsg = 'ERROR EXIT CODE'
    commands.getoutput('cat '+OUTFILE+' >> '+LOGFILE)    # append to LOGFILE
    commands.getoutput('cat '+OUTFILE+' >> '+ERRORFILE)  # append to ERRORFILE
  else:
    res = check_residual(case,residual)
    if(res[0] == 1):
      smsg = 'OK  '
      commands.getoutput('cat '+OUTFILE+' >> '+LOGFILE)    # append to LOGFILE
    else:
      smsg = 'FAIL'
      lmsg = 'Expected '+str(residual)+' Got '+str(res[1])
      commands.getoutput('cat '+OUTFILE+' >> '+LOGFILE)    # append to LOGFILE
      commands.getoutput('cat '+OUTFILE+' >> '+ERRORFILE)  # append to ERRORFILE

  print_stdout(smsg+" : "+time+" : "+case+" : "+lmsg+"\n")

# --- function to check the residual was correct ------
def check_residual(case,trueres):
  convfile = open(convname, 'r')               # open the file
  contents = convfile.read()                   # get the contents of the convergence file
  convfile.close()                             # close the file

  lines = contents.splitlines()                # get the lines in the convergence file
  lastline = lines[len(lines)-1]               # get the last line
  fields = string.splitfields(lastline)        # separate the fields in the last line

  calcres = string.atof(fields[1])             # get the calculted residual

  if (abs(calcres - trueres) < abs(fuzz*calcres)):
    return (1,calcres)
  else:
    return (0,calcres)

# --- function to run the testcase -------------------
def run_testcase(case,residual):
  stp = StopWatch()
  stp.Start()
  output = commands.getstatusoutput(PROG + " " + case + "&>" + OUTFILE)
  stp.Stop()

  print_output(case,stp,output,residual)

# --- function to run the unit tests -------------------
def run_unittest(unittest):
  output = commands.getstatusoutput("$COOLFLUID_BUILD_DIR/" + unittest + " &> " + OUTFILE)

  if(output[0] != 0):
        print_stdout(unittest+' : FAIL\n')
  else:
        print_stdout(unittest+' : OK\n')

  commands.getoutput('cat '+OUTFILE+' >> '+LOGFILE)    # append to LOGFILE

# --- function to print the enviromental variables -------------------
def printenv():
  s = ''
  for e in os.environ.items():
    s = s + e[0] + '=' + repr(e[1]) + ' '
  return s

# --- main body  -------------------------------------

output = commands.getoutput("rm -f " + LOGFILE)
output = commands.getoutput("rm -f " + ERRORFILE)

testcases = open(CATALOGFILE, 'r')           # open the testcases catalog file
testcont = testcases.read()                  # get the contents of the convergence file
testcases.close()                            # close the file

testlines = testcont.splitlines()            # get the lines in the testcase file

t = datetime.datetime.now()

# prints the relevant environmental variables
print_stdout("---=== Checking COOLFluiD Testcases ===--- \n\n")
print_stdout("Date: "+t.strftime("%A, %d %b %Y @ %Hh%M")+"\n\n")
print_stdout("Environmental variables\n")
commands.getstatusoutput("echo $COOLFLUID_DIR")
print_stdout("Current dir : "+os.environ['PWD']+"\n")
print_stdout("COOLFLUID_DIR : "+os.environ['COOLFLUID_DIR']+"\n")
print_stdout("COOLFLUID_BUILD_DIR : "+os.environ['COOLFLUID_BUILD_DIR']+"\n")
print_stdout("COOLFLUID_PROG : "+os.environ['COOLFLUID_PROG']+"\n")

globalstp = StopWatch()
globalstp.Start()

# run the testcases and compare residuals
print_stdout("\nRunning testcases\n")

for line in testlines:
  fields  = string.splitfields(line)         # separate the fields in the last line
  case    = fields[0]                        # get the case filename
  trueres = string.atof(fields[1])           # get the expected residual
  run_testcase(case,trueres)

# run the unit tests
print_stdout("\nRunning unit tests\n")
run_unittest("test/math/testMath")
run_unittest("test/shapefunctions/testShapeFunctions")
run_unittest("test/utils/testUtils")
run_unittest("test/common/testCommon")
run_unittest("test/framework/testFramework")

globalstp.Stop()

totaltime = globalstp.getStr()
print_stdout("\nTotal time: "+totaltime+"\n")

print "\nOutput written to file: "+CHECKLOG

# clean up files
output = commands.getoutput("rm -f " + OUTFILE)
