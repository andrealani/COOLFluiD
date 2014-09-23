#ifndef COOLFluiD_cient_GlobalLog_h
#define COOLFluiD_cient_GlobalLog_h

//////////////////////////////////////////////////////////////////////////////

#include <QObject>

class QString;

namespace COOLFluiD
{
  namespace client
  {
    
    //////////////////////////////////////////////////////////////////////////////
    
    /// @brief Allows to write to a log, from everywhere.
    
    /// This class is a singleton. It provides static and non-static methods
    /// the write a message or an error message. Static methods are provided for
    /// convenience. Doing 
    /// @code GlobalLog::message("Hello world!"); @endcode is equivalent to
    /// @code GlobalLog::getInstance()->writeMessage("Hello world!"); @endcode
    /// Each time an error or a message is written, the appropriate signal
    /// is emitted with the message.
    
    /// @author Quentin Gasper
    
    class GlobalLog : public QObject
    {
      Q_OBJECT
      
    public:
      
      /// @brief Gives a pointer to the class unique m_instance.
      
      /// If no m_instance exists, a new one is created.
      /// @return Returns a pointeur to the class unique m_instance.
      static GlobalLog * getInstance();
      
      /// @brief Writes an error message.
      
      /// If @c msg is empty, nothing is done.
      /// @param msg The message to write.
      /// @param fromServer If @c true, the message comes from the server.
      void writeError(const QString & msg, bool fromServer = false) const;
      
      /// @brief Writes a message.
      
      /// If @c msg is empty, nothing is done.
      /// @param msg The message to write.
      /// @param fromServer If @c true, the message comes from the server.
      void writeMessage(const QString & msg, bool fromServer = false) const;
      
      /// @brief Writes an error message.
      
      /// This is an overloaded function, provided for convenience.
      /// If @c msg is empty, nothing is done.
      /// @param msg The message to write.
      /// @param fromServer If @c true, the message comes from the server.
      /// @see writeError
      static void error(const QString & msg, bool fromServer = false);
      
      /// @brief Writes a message.
      
      /// This is an overloaded function, provided for convenience.
      /// If @c msg is empty, nothing is done.
      /// @param msg The message to write.
      /// @param fromServer If @c true, the message comes from the server.
      /// @see writeMessage
      static void message(const QString & msg, bool fromServer = false);
      
    signals:
      
      /// @brief Signal emitted when an errror message has to be write in the log.
      
      /// @param msg The error message
      /// @param fromServer If @c true, the message comes from the server.
      void sigError(const QString & msg, bool fromServer) const;
      
      /// @brief Signal emitted when a message has to be write in the log.
      
      /// @param msg The message
      /// @param fromServer If @c true, the message comes from the server.
      void sigMessage(const QString & msg, bool fromServer) const;
      
    private:
      
      /// @brief Class unique m_instance
      static GlobalLog * m_instance;
      
      /// @brief Constructor
      GlobalLog();
      
      
    }; // class GlobalLog
    
    //////////////////////////////////////////////////////////////////////////////
    
  } // namespace client
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif
