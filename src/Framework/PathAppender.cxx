// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/PE.hh"
#include "Framework/SubSystemStatus.hh"
#include "Framework/SimulationStatus.hh"
#include "Framework/PathAppender.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;
using namespace COOLFluiD::Common;

namespace COOLFluiD {

  namespace Framework {

//////////////////////////////////////////////////////////////////////////////

PathAppender::PathAppender() :
  m_appendIterFlag(false),
  m_appendTimeFlag(false)
{
  m_appendParallelFlag = PE::GetPE().IsParallel();
}

//////////////////////////////////////////////////////////////////////////////

PathAppender::~PathAppender()
{
}

//////////////////////////////////////////////////////////////////////////////

PathAppender& PathAppender::getInstance()
{
  static PathAppender appender;
  return appender;
}

//////////////////////////////////////////////////////////////////////////////

boost::filesystem::path
PathAppender::appendParallel(const boost::filesystem::path& fpath) const
{
  using namespace boost::filesystem;

  ostringstream add;
  if (m_appendParallelFlag) {
    add << "-P" << PE::GetPE().GetRank("Default");
  }

  return fpath.branch_path() / ( basename(fpath) + add.str() + extension(fpath) );
}

//////////////////////////////////////////////////////////////////////////////

boost::filesystem::path
PathAppender::appendIter(const boost::filesystem::path& fpath) const
{
  using namespace boost::filesystem;

  ostringstream add;
  if (m_appendIterFlag) {
    CFuint iter = SubSystemStatusStack::getActive()->getNbIter();
    add << "-iter_" << iter;

    if(SubSystemStatusStack::getActive()->doingSubIterations()) add << "-sub_" << SubSystemStatusStack::getActive()->getSubIter();
  }

  return fpath.branch_path() / ( basename(fpath) + add.str() + extension(fpath) );
}

//////////////////////////////////////////////////////////////////////////////

boost::filesystem::path
PathAppender::appendGlobalIter(const boost::filesystem::path& fpath) const
{
  using namespace boost::filesystem;

  ostringstream add;
  const bool appendGlobalIter = SimulationStatus::getInstance().isAppendIter();
  if (appendGlobalIter) {
    CFreal iter = SimulationStatus::getInstance().getNbIter();
    add << "-globaliter_" << iter;
  }

  return fpath.branch_path() / ( basename(fpath) + add.str() + extension(fpath) );
}

//////////////////////////////////////////////////////////////////////////////

boost::filesystem::path
PathAppender::appendTime(const boost::filesystem::path& fpath) const
{
  using namespace boost::filesystem;

  ostringstream add;
  if (m_appendTimeFlag) {

    CFreal time = SubSystemStatusStack::getActive()->getCurrentTimeDim();

    if(SubSystemStatusStack::getActive()->doingSubIterations()
     && !SubSystemStatusStack::getActive()->isSubIterationLastStep() )
    {
        time += SubSystemStatusStack::getActive()->getDTDim();
    }

    std::string time_string = Common::StringOps::to_str(time);
    Common::StringOps::subst(".", "_",time_string);
    add << "-time_" << time_string ;
  }

  return fpath.branch_path() / ( basename(fpath) + add.str() + extension(fpath) );
}

//////////////////////////////////////////////////////////////////////////////

boost::filesystem::path
PathAppender::appendAllInfo(const boost::filesystem::path& fpath, 
                            bool appendIterFlag, bool appendTimeFlag,
                            bool appendParallelFlag)
{
  m_appendIterFlag = appendIterFlag;
  m_appendTimeFlag = appendTimeFlag;
  m_appendParallelFlag = appendParallelFlag;

  boost::filesystem::path newPath =
    appendTime(
    appendGlobalIter(
    appendIter(
    appendParallel(fpath) ) ) ) ;

  //set back to default
  m_appendIterFlag = false;
  m_appendTimeFlag = false;
  m_appendParallelFlag = PE::GetPE().IsParallel();

  return newPath;
}

//////////////////////////////////////////////////////////////////////////////

boost::filesystem::path
PathAppender::appendAllInfo(const boost::filesystem::path& fpath)
{
  m_appendIterFlag = true;
  m_appendTimeFlag = true;

  boost::filesystem::path newPath =
    appendTime(
    appendGlobalIter(
    appendIter(
    appendParallel(fpath) ) ) ) ;

  //set back to default
  m_appendIterFlag = false;
  m_appendTimeFlag = false;

  return newPath;
}

//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

boost::filesystem::path
PathAppender::appendCustom(const boost::filesystem::path& fpath, std::string& custom) const
{
  using namespace boost::filesystem;

  ostringstream add;
  add << "-" << custom;

  return fpath.branch_path() / ( basename(fpath) + add.str() + extension(fpath) );
}

//////////////////////////////////////////////////////////////////////////////

  } // namespace Framework

} // namespace COOLFluiD
