// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <fstream>

#include "Common/CFLog.hh"
#include "Common/StringOps.hh"
#include "Common/FilesystemException.hh"
#include "Config/ConfigArgs.hh"
#include "Config/ConfigFileReader.hh"
#include "Config/NestedConfigFileReader.hh"
#include "Config/BadMatchException.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;

namespace COOLFluiD {

  namespace Config {

//////////////////////////////////////////////////////////////////////////////

const std::string NestedConfigFileReader::SEPARATOR   = "=";
const CFchar   NestedConfigFileReader::CONTINUATOR = '\\';
const std::string NestedConfigFileReader::COMMENTOR   = "#";
const std::string NestedConfigFileReader::METACOMMENTOR   = "###";
const std::string NestedConfigFileReader::BLOCKCOMMENTORBEGIN   = "!#";
const std::string NestedConfigFileReader::BLOCKCOMMENTOREND   = "#!";

//////////////////////////////////////////////////////////////////////////////

NestedConfigFileReader::NestedConfigFileReader ()
{
}

//////////////////////////////////////////////////////////////////////////////

NestedConfigFileReader::~NestedConfigFileReader()
{
}

//////////////////////////////////////////////////////////////////////////////

void NestedConfigFileReader::parse (const std::string& filename_in, ConfigArgs& args)
{
  if(!filename_in.empty()) {

  std::string filename = expandFileName(filename_in);

  std::ifstream iss(filename.c_str());
  if(!iss) throw Common::FilesystemException (FromHere(),"Could not open file: " + filename);

  bool commentedBlock = false;
  while (iss)
  {

    std::string line;
    getline(iss, line);

    // Commented blocks are removed from the file reading
    if(commentedBlock == true) {
      std::string::size_type blockendcmtpos = line.find(BLOCKCOMMENTOREND);
      if (blockendcmtpos != std::string::npos) {
        // strip to the end of the line
        commentedBlock = false;
        line.erase(0,blockendcmtpos + BLOCKCOMMENTOREND.size());
      }
      else {
        line.erase();
      }
    }

    // Check for meta-comment with include directive
    std::string::size_type metacmtpos = line.find(METACOMMENTOR);
    if ( metacmtpos != std::string::npos )
    {
      std::vector<std::string> wds = Common::StringOps::getWords(line);

      if (wds.size() > 1 ) // more than metacomment is present
      {
         // check for include directive
          if ( wds[1] == "IncludeCase" )
          {
            if ( wds.size() > 2 )
            {
              ConfigFileReader other_file_reader;
              other_file_reader.parse( wds[2], args);
            }
            else
              CFLogWarn ( "IncludeCase directive has no file defined\n" );
          }
      }
    }

    // Check for start of a commented block
    std::string::size_type blockbegincmtpos = line.find(BLOCKCOMMENTORBEGIN);
    if (blockbegincmtpos != std::string::npos) {
      // strip to the end of the line
      line.erase(blockbegincmtpos);
      commentedBlock = true;
    }

    /// @warning Weak comment stripper. It fails on comment characters inside a std::string.
    std::string::size_type cmtpos = line.find(COMMENTOR);
    if (cmtpos != std::string::npos) {
      // strip to the end of the line
      line.erase(cmtpos);
    }

    // strip trailing whitespaces
    Common::StringOps::trimRear(line);

  if(!line.empty())
  {

//     std::string::size_type cont = line.rfind(CONTINUATOR);
//     if(cont != std::string::npos) {
//       line = std::string(line,0,cont);
//       std::string nextline;
//       getline(iss, nextline);
//       nextline.trimRear();
//       line += nextline;
//     }

    // if the last character is a backslash, keep concatenating the lines
    while (iss && line[line.length() - 1] == CONTINUATOR) {
      line = std::string(line, 0, line.length() - 1);
      std::string nextline;
      getline(iss, nextline);
      Common::StringOps::trimRear(nextline);
      line += nextline;
    }

    CFint shift = 0;
    size_t eqpos = line.find(SEPARATOR);
    if (eqpos != std::string::npos) {
      std::string preSeparator(line, eqpos-1, 1);
      if(preSeparator =="+"){
        shift = 1;
      }
      std::string label(line, 0, eqpos-shift);
      std::string value(line, eqpos + 1);

      Common::StringOps::trim(label); // trim front and back around the label
      Common::StringOps::trim(value); // they might have leading and tailing whitespace

      if((preSeparator =="+") && (args[label] != ""))
      {
        args[label] = args[label] + " " + value;
      }
      else
      {
        if(args[label] != "") throw Config::BadMatchException (FromHere(),"Overriding previously specified option for: " + label);
        args[label] = value;
      }
    }
  } // not empty
  } // while loop
  } // filename is empty
}

//////////////////////////////////////////////////////////////////////////////

std::string NestedConfigFileReader::expandFileName(const std::string& fname) const
{
  std::string exp = fname;
  if (fname[0] == '~') {
    char* home = getenv("HOME");
    if (home) {
      exp = std::string(home) + fname.substr(1);
    }
  }
  return exp;
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Config

} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////
