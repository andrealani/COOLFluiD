#include "logcpp/PortabilityImpl.hh"

#include "ClientServer/server/RemoteClientAppender.hh"

namespace COOLFluiD {
  namespace server {
    
    RemoteClientAppender::RemoteClientAppender(const std::string& name) : logcpp::LayoutAppender(name)
    {}
    
    RemoteClientAppender::~RemoteClientAppender()
    {
      close();
    }
    
    void RemoteClientAppender::close()
    {
      // empty
    }
    
    void RemoteClientAppender::_append(const logcpp::LoggingEvent& event)
    {
      emit newData( _getLayout().format(event).c_str() );
    }
    
    bool RemoteClientAppender::reopen()
    {
      return true;
    }
    
  } // server
} // coolfluid
