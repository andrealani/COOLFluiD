// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include <cctype> // for to_lower and to_upper

#include "Common/StringOps.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std; // for all the std::string stuff

namespace COOLFluiD {
namespace Common {

//////////////////////////////////////////////////////////////////////////////

void StringOps::join (const std::string& delim, std::string parts[], unsigned int nparts, std::string& out)
{
  out = parts[0];
  for (size_t i = 1; i < nparts; ++i)
  {
    out += delim;
    out += parts[i];
  }
}

//////////////////////////////////////////////////////////////////////////////

void StringOps::toLower (std::string& out)
{
  const string::iterator end = out.end();
  for (string::iterator itr = out.begin(); itr != end; ++itr) {
    *itr = tolower(*itr);
  }
}

//////////////////////////////////////////////////////////////////////////////

void StringOps::toUpper (std::string& out)
{
  const string::iterator end = out.end();
  for (string::iterator itr = out.begin(); itr != end; ++itr) {
    *itr = toupper(*itr);
  }
}

//////////////////////////////////////////////////////////////////////////////

void StringOps::subst (const std::string& lhs, const std::string& rhs, std::string& out)
{
  std::string::size_type i = 0;
  while (i < out.length())
  {
    std::string substr (out, i);
    std::string::size_type idx = out.find(lhs, i);
    if (idx == std::string::npos)
    {
      // we're at the end
      break;
    }
    else
    {
      out.replace(idx, lhs.length(), rhs);
      i = idx + rhs.length();
      if (!rhs.length())
         ++i;
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void StringOps::trimFront (std::string& out)
{
  size_t startpos = out.find_first_not_of(" \t");
  if( string::npos != startpos )
        out = out.substr( startpos );
}

//////////////////////////////////////////////////////////////////////////////

void StringOps::trimRear (std::string& out)
{
  size_t endpos = out.find_last_not_of(" \t");
  if( string::npos != endpos )
        out = out.substr( 0, endpos+1 );
}

//////////////////////////////////////////////////////////////////////////////

void StringOps::trim (std::string& out)
{
  size_t startpos = out.find_first_not_of(" \t");
  size_t endpos = out.find_last_not_of(" \t");
  if(( string::npos == startpos ) || ( string::npos == endpos))
  {  out = ""; }
  else
     out = out.substr( startpos, endpos-startpos+1 );
}

//////////////////////////////////////////////////////////////////////////////

void StringOps::trim2 (std::string& out)
{
  size_t startpos = out.find_first_not_of(" \t\n");
  size_t endpos = out.find_last_not_of(" \t\n");
  if(( string::npos == startpos ) || ( string::npos == endpos))
    {  out = ""; }
  else
     out = out.substr( startpos, endpos-startpos+1 );
}

//////////////////////////////////////////////////////////////////////////////

std::vector<std::string> StringOps::getWords (const std::string& in)
{
  std::string s = in;            // copy
  vector<std::string> words;     // return vector

  bool inWord = false;           // whether we're in a word
  std::string word;              // current word
  for (size_t i = 0, len = s.length(); i < len; ++i)
  {
    CFchar ch = s[i];
    if (inWord)
    {
      if (isspace(ch))
      {
        words.push_back(word);
        word = "";
        inWord = false;
      }
      else
      {
        word += ch;
      }
    }
    else
      if (!isspace(ch))
      {
        word += ch;
        inWord = true;
      }
  }

  // grab the last one, if any.
  if (inWord)
    words.push_back(word);

  return words;
}

//////////////////////////////////////////////////////////////////////////////

vector<std::string> StringOps::getWords  (const std::string& in, const CFchar sep )
{
  std::string s = in;                // copy
  vector<std::string> words;         // return vector

  bool inWord = false;        // whether we're in a word
  std::string word;                // current word
  for (size_t i = 0, len = s.length(); i < len; ++i) {
    CFchar ch = s[i];
    if (inWord) {
      if (ch == sep) {
        words.push_back(word);
        word = "";
        inWord = false;
      }
      else {
        word += ch;
      }
    }
    else if (ch != sep) {
      word += ch;
      inWord = true;
    }
  }

  // grab the last one, if any.
  if (inWord)
    words.push_back(word);

  return words;
}

//////////////////////////////////////////////////////////////////////////////

bool StringOps::startsWith (const std::string& in, const std::string& str)
{
  return in.length() >= str.length() && in.substr(0, str.length()) == str;
}

//////////////////////////////////////////////////////////////////////////////

bool StringOps::endsWith (const std::string& in, const std::string& str)
{
 return in.length() >= str.length() && in.substr(in.length() - str.length()) == str;
}

//////////////////////////////////////////////////////////////////////////////

} // namespace Common
} // namespace COOLFluiD
