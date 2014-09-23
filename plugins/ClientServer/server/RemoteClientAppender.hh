#ifndef COOLFluiD_server_RemoteClientAppender_hh
#define COOLFluiD_server_RemoteClientAppender_hh

#include <string>
#include <iostream>

#include <QObject>

#include <logcpp/Portability.hh>
#include <logcpp/LayoutAppender.hh>

namespace COOLFluiD {
  namespace server {
    
    /// Appends LoggingEvents to the remote client log window.
    
    class RemoteClientAppender : public QObject, public logcpp::LayoutAppender
    {
      Q_OBJECT
      
    public:
      
      RemoteClientAppender(const std::string& name);
      virtual ~RemoteClientAppender();
      
      virtual bool reopen();
      virtual void close();
      
    protected:
      virtual void _append(const logcpp::LoggingEvent& event);
      
    signals:
      void newData(const QString & m_data);
    };
    
  } // server
} // coolfluid

#endif // COOLFluiD_server_RemoteClientAppender_hh
