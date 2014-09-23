/*
 * SimpleConfigurator.cpp
 *
 * Copyright 2001, Glen Scott. All rights reserved.
 *
 * See the COPYING file for the terms of usage and distribution.
 */
#include "PortabilityImpl.hh"

#ifdef CF_HAVE_UNISTD_H
  #include <unistd.h>
#endif
#ifdef LOGCPP_HAVE_IO_H
#    include <io.h>
#endif

#include <stdio.h>
#include <iostream>
#include <string>
#include <fstream>
#include <stdio.h>

#include <logcpp/Category.hh>
#include <logcpp/Appender.hh>
#include <logcpp/OstreamAppender.hh>
#include <logcpp/FileAppender.hh>
#include <logcpp/RollingFileAppender.hh>
#include <logcpp/Layout.hh>
#include <logcpp/BasicLayout.hh>
#include <logcpp/SimpleLayout.hh>
#include <logcpp/Priority.hh>
#include <logcpp/NDC.hh>
#include <logcpp/PatternLayout.hh>
#include <logcpp/SimpleConfigurator.hh>
// #include <logcpp/RemoteSyslogAppender.hh>

namespace logcpp {

    void SimpleConfigurator::configure(const std::string& initFileName) throw (ConfigureFailure) {
        std::ifstream initFile(initFileName.c_str());

        if (!initFile) {
            throw ConfigureFailure(std::string("Config File ") + initFileName + " does not exist or is unreadable");
        }

        configure(initFile);
    }

    void SimpleConfigurator::configure(std::istream& initFile) throw (ConfigureFailure) {
        std::string nextCommand;
        std::string categoryName;

        while (initFile >> nextCommand) {
            /* skip comment lines */
            if (nextCommand[0] == '#') {
                std::string dummy;
                std::getline(initFile, dummy);
                continue;
            }
            /* stop on missing categoryName */
            if (!(initFile >> categoryName))
                break;

            logcpp::Category& category =
                (categoryName.compare("root") == 0) ?
                logcpp::Category::getRoot() :
                logcpp::Category::getInstance(categoryName);

            if (nextCommand.compare("appender") == 0) {
                std::string layout;
                std::string appenderName;

                if (initFile >> layout >> appenderName) {
                    logcpp::Appender* appender;
                    if (appenderName.compare("file") == 0) {
                        std::string logFileName;
                        if (!(initFile >> logFileName)) {
                            throw ConfigureFailure("Missing filename for log file logging configuration file for category: " + categoryName);
                        }
                        appender = new logcpp::FileAppender(categoryName, logFileName);
                    }
                    else if (appenderName.compare("rolling") == 0) {
                        std::string logFileName;
                        size_t maxFileSize;
                        unsigned int maxBackupIndex=1;
                        if (!(initFile >> logFileName)) {
                            throw ConfigureFailure("Missing filename for log file logging configuration file for category: " + categoryName);
                        }
				if (!(initFile >> maxFileSize)) {
                            throw ConfigureFailure("Missing maximum size for log file logging configuration file for category: " + categoryName);
                        }
                        if (!(initFile >> maxBackupIndex)) {
                            throw ConfigureFailure("Missing maximum backup index for log file logging configuration file for category: " + categoryName);
                        }
                        appender = new logcpp::RollingFileAppender(categoryName, logFileName, maxFileSize, maxBackupIndex);
                    }
                    else if (appenderName.compare("console") == 0) {
                        appender =
                            new logcpp::OstreamAppender(categoryName, &std::cout);
                    }
                    else if (appenderName.compare("stdout") == 0) {
                        appender =
                            new logcpp::FileAppender(categoryName, ::dup(fileno(stdout)));
                    }
                    else if (appenderName.compare("stderr") == 0) {
                        appender =
                            new logcpp::FileAppender(categoryName, ::dup(fileno(stderr)));
                    }
#if 0 // no RemoteSysLogAppender ( must fix build in Win32 )
                    else if (appenderName.compare("remotesyslog") == 0) {
                        std::string syslogName;
                        std::string relayer;
                        int facility;
                        int portNumber;
                        if (!(initFile >> syslogName)) {
                            throw ConfigureFailure("Missing syslogname for SysLogAppender for category: " + categoryName);
                        }
                        if (!(initFile >> relayer)) {
                            throw ConfigureFailure("Missing syslog host for SysLogAppender for category: " + categoryName);
                        }
                        if (!(initFile >> facility)) {
                            facility = LOG_USER;
                        }
                        if (!(initFile >> portNumber)) {
                            portNumber = 514;
                        }
                        appender =
                            new logcpp::RemoteSyslogAppender(categoryName, syslogName, relayer, facility, portNumber);
                    }
#endif
                    else {
                        throw ConfigureFailure("Invalid appender name (" +
                                               appenderName +
                                               ") in logging configuration file for category: " +
                                               categoryName);
                    }
                    if (layout.compare("basic") == 0)
                        appender->setLayout(new logcpp::BasicLayout());
                    else if (layout.compare("simple") == 0)
                        appender->setLayout(new logcpp::SimpleLayout());
                    else if (layout.compare("pattern") == 0) {
                        logcpp::PatternLayout *layout =
                            new logcpp::PatternLayout();
			initFile >> std::ws; // skip whitespace
                        char pattern[1000];
                        initFile.getline(pattern, 1000);
                        layout->setConversionPattern(std::string(pattern));
			appender->setLayout(layout);
                    }
                    else {
                        throw ConfigureFailure("Invalid layout (" + layout +
                                               ") in logging configuration file for category: " +
                                               categoryName);
                    }
                    category.addAppender(appender);
                }
            }
            else if (nextCommand.compare("priority") == 0) {
                std::string priority;
                if (!(initFile >> priority)) {
                    throw ConfigureFailure("Missing priority in logging configuration file for category: " + categoryName);
                }

                try {
                    category.setPriority(logcpp::Priority::getPriorityValue(priority));
                } catch(std::invalid_argument) {
                    throw ConfigureFailure("Invalid priority ("+priority+") in logging configuration file for category: "+categoryName);
                }
            }
            else if (nextCommand.compare("category") == 0) {
                /*
                  This command means we should "refer" to the category
                  (in order to have it created). We've already done this
                  in common setup code for all commands.
                */
            }
            else {
                throw ConfigureFailure("Invalid format in logging configuration file. Command: " + nextCommand);
            }
        }
    }
}



