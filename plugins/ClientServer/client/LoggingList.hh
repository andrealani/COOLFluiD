#ifndef COOLFluiD_client_LoggingList_h
#define COOLFluiD_client_LoggingList_h

#include <QTextEdit>
#include <QHash>

class QColor;

/////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD
{
  namespace client
  {
    
    /////////////////////////////////////////////////////////////////////////////
    
    /// @brief Defines available message types for the log.
    enum LogMessageType
    {
      /// @brief The associated message is a normal message from the client.
      TYPE_NORMAL,
      /// @brief The associated message is a error message from the client.
      TYPE_ERROR,
      /// @brief The associated message is a error message from the server.
      TYPE_SERVER_ERROR,
      /// @brief The associated message is a normal message from the server.
      TYPE_SERVER_NOTICE,
      /// @brief The associated message is a informative message from the server.
      TYPE_SERVER_INFO,
      /// @brief The associated message is a warning message from the server.
      TYPE_SERVER_WARN
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
    /// @brief Manages a graphical log component.
    
    class LoggingList : public QTextEdit
    {
      Q_OBJECT 
      
    public:
      
      /// @brief Constructor
      
      /// @param parent Parent. May be @c NULL.
      /// @param maxLogLines Number of lines before the log must be cleared. If 0,
      /// the log is never cleared.
      LoggingList(QWidget * parent = NULL, unsigned int maxLogLines = 100000);
      
      /// @brief Destructor
      
      /// Frees all allocated memory. Parent is not destroyed.
      ~LoggingList();
      
      /// @brief Assigns a new color to a message type.
      
      /// @param type Message type.
      /// @param color Color. It must be a valid color (@c color.isValid() must 
      /// return @c true ).
      /// @throw IllegalValueException If the color is not valid.
      void setColor(LogMessageType type, const QColor & color);
      
      /// @brief Gives a color associated to a type.
      
      /// @param type Type of which the color is wanted.
      /// @return Returns the color.
      QColor getColor(LogMessageType type) const;
      
      /// @brief Sets the maximum number of lines before the log is cleared.
      
      /// @param maxLogLines Maximum number of lines. If 0, the log is never 
      /// cleared.
      /// @warning Be careful when modifying this value. The more lines the log
      /// contains, the more it takes memory to store its content.
      void setMaxLogLines(unsigned int maxLogLines = 100000);
      
      /// @brief Gives the maximum number of lines before the log is cleared.
      
      /// @return Returns the maximum number of lines before the log is cleared.
      unsigned int getMaxLogLines() const;
      
      /// @brief Gives the number of lines the log contains.
      
      /// @return Returns the number of lines the log contains.
      unsigned int getLogLinesCounter() const;
      
      /// @brief Appends a message to the log.
      
      /// To save memory, the list is cleared every @c #maxLogLines lines.
      /// @param string Message to append.
      /// @param type Type of message. Must be one the values defined by 
      /// @c LogMessageType enum.    
      void appendToLog(const QString & string, LogMessageType type);
      
      public slots:
      
      void message(const QString & message, bool fromServer);
      
      void error(const QString & message, bool fromServer);
      
      /// @brief Clears the log
      
      /// Use this method instead of the base class method @c clear() to ensure
      /// that the line counter is reset.
      void clearLog();
      
      
    private:
      
      /// @brief Colors used for different message types.
      
      /// The key is the type, the value is the color.
      QHash<LogMessageType, QColor> m_logColors;
      
      /// @brief Prefixes prepended to message depending on their type.
      
      /// The key is the type, the value is the prefix.
      QHash<LogMessageType, QString> m_prefixes;
      
      /// @brief Number of lines in the log.
      unsigned int m_logLinesCounter;
      
      /// @brief Number of lines before the log must be cleared.
      
      /// If 0, the log is never cleared.
      unsigned int m_maxLogLines;
      
      
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
  } // namespace client
} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_client_LoggingList_h
