LIST(APPEND logcpp_files
PropertyConfigurator.hh
BasicConfigurator.hh
Properties.hh
NDC.hh
StringQueueAppender.hh
Threading.hh
AppenderSkeleton.hh
RollingFileAppender.hh
SimpleConfigurator.hh
CategoryStream.hh
Configurator.hh
BoostThreads.hh
BoostThreads.cpp
TimeStamp.hh
Layout.hh
FixedContextCategory.hh
StringUtil.hh
Portability.hh
PatternLayout.hh
LoggingEvent.hh
Priority.hh
Appender.hh
# RemoteSyslogAppender.hh
Category.hh
HierarchyMaintainer.hh
PropertyConfiguratorImpl.hh
LayoutAppender.hh
FileAppender.hh
SimpleLayout.hh
Filter.hh
BasicLayout.hh
Export.hh
OstreamAppender.hh
# DummyThreads.hh
PortabilityImpl.hh
HierarchyMaintainer.cpp
PatternLayout.cpp
OstreamAppender.cpp
PropertyConfigurator.cpp
# RemoteSyslogAppender.cpp
Configurator.cpp
NDC.cpp
# DummyThreads.cpp
LoggingEvent.cpp
Category.cpp
Priority.cpp
StringQueueAppender.cpp
Filter.cpp
RollingFileAppender.cpp
TimeStamp.cpp
SimpleConfigurator.cpp
Appender.cpp
PortabilityImpl.cpp
FileAppender.cpp
BasicConfigurator.cpp
PropertyConfiguratorImpl.cpp
Properties.cpp
BasicLayout.cpp
AppenderSkeleton.cpp
FixedContextCategory.cpp
SimpleLayout.cpp
LayoutAppender.cpp
CategoryStream.cpp
StringUtil.cpp
config.h
config-win32.h
)

LIST ( APPEND OPTIONAL_dirfiles RemoteSyslogAppender.hh RemoteSyslogAppender.cpp )

LIST ( APPEND logcpp_libs   ${CF_Boost_LIBRARIES} )

CF_ADD_KERNEL_LIBRARY ( logcpp )

CF_WARN_ORPHAN_FILES()
