// Copyright (C) 2012 von Karman Institute for Fluid Dynamics, Belgium
//
// This software is distributed under the terms of the
// GNU Lesser General Public License version 3 (LGPLv3).
// See doc/lgpl.txt and doc/gpl.txt for the license text.

#include "Common/CFLog.hh"

#include "logcpp/FileAppender.hh"
#include "logcpp/OstreamAppender.hh"
#include "logcpp/PatternLayout.hh"

//////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace COOLFluiD::Common;

namespace COOLFluiD {

//////////////////////////////////////////////////////////////////////////////

CFLogger& CFLogger::getInstance ()
{
  static CFLogger logger;
  return logger;
}

//////////////////////////////////////////////////////////////////////////////

CFLogger::CFLogger ()
{
  // std::cout << " ++++++ Constructing CFLogger ++++++ " << std::endl;

  // set basic output policies
  logcpp::Category::getRoot().removeAllAppenders();
  getTraceLogger().removeAllAppenders();
  getMainLogger().removeAllAppenders();

  logcpp::PatternLayout* stdout_layout = new logcpp::PatternLayout();
  std::string stdout_format("%m");
  stdout_layout->setConversionPattern(stdout_format);

  logcpp::Appender* stdout_appender = new logcpp::OstreamAppender("StdoutAppender", &std::cout);
  stdout_appender->setLayout(stdout_layout);

  getMainLogger().addAppender(stdout_appender);
}

//////////////////////////////////////////////////////////////////////////////

CFLogger::~CFLogger ()
{
  // std::cout << " ++++++ Destroying CFLogger ++++++ " << std::endl;
}

//////////////////////////////////////////////////////////////////////////////

logcpp::Category& CFLogger::getMainLogger()
{
  static logcpp::Category& logger = logcpp::Category::getInstance("COOLFluiD Logger");
  return logger;
}

//////////////////////////////////////////////////////////////////////////////

logcpp::Category& CFLogger::getTraceLogger()
{
  static logcpp::Category& tracer = logcpp::Category::getInstance("COOLFluiD Tracer");
  return tracer;
}

//////////////////////////////////////////////////////////////////////////////

void CFLogger::setMainLoggerLevel(CFuint level)
{
  CFuint mylevel = level;
  if(level > DEBUG_MAX) mylevel = DEBUG_MAX;
  if(level < ERROR)     mylevel = ERROR;
  getMainLogger().setPriority(static_cast<CFLogLevel>(mylevel));
}

//////////////////////////////////////////////////////////////////////////////

void CFLogger::setTraceToStdOut(bool traceToStdOut)
{
  if (traceToStdOut)
  {
    logcpp::PatternLayout* stdout_layout = new logcpp::PatternLayout();
    std::string stdout_format("%m");
    stdout_layout->setConversionPattern(stdout_format);

    logcpp::Appender* ts_appender = new logcpp::OstreamAppender("StdoutAppender", &std::cout);
    ts_appender->setLayout(stdout_layout);

    getTraceLogger().addAppender(ts_appender);
  }
}

//////////////////////////////////////////////////////////////////////////////

} // namespace COOLFluiD
