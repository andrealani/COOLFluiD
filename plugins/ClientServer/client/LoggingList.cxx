#include <QtGui>

#include "ClientServer/client/LoggingList.hh"

using namespace COOLFluiD::client;

LoggingList::LoggingList(QWidget * parent, unsigned int maxLogLines)
: QTextEdit(parent)
{
  m_maxLogLines = maxLogLines;
  
  this->setReadOnly(true);
  this->setWordWrapMode(QTextOption::NoWrap);
  
  m_prefixes[ TYPE_SERVER_NOTICE ] = "[ SERVER &nbsp;] ";
  m_prefixes[ TYPE_SERVER_INFO ] =   "[ INFO &nbsp; &nbsp;] ";
  m_prefixes[ TYPE_SERVER_WARN ] =   "[ WARNING ] ";
  m_prefixes[ TYPE_SERVER_ERROR ] =  "[ ERROR &nbsp; ] ";
  m_prefixes[ TYPE_ERROR ]        = "[ CLIENT &nbsp;]";
  m_prefixes[ TYPE_NORMAL ]       =  "[ CLIENT &nbsp;]";
  
  m_logColors[ TYPE_SERVER_NOTICE ] = QColor("black");
  m_logColors[ TYPE_SERVER_INFO ] =   QColor("black");
  m_logColors[ TYPE_SERVER_WARN ] =   QColor("darkorange");
  m_logColors[ TYPE_SERVER_ERROR ] =  QColor("red");
  m_logColors[ TYPE_ERROR ] =         QColor("red");
  m_logColors[ TYPE_NORMAL ] =        QColor("black");
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

LoggingList::~LoggingList()
{
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void LoggingList::setColor(LogMessageType type, const QColor & color)
{
  if(color.isValid())
    m_logColors[type] = color.name();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

QColor LoggingList::getColor(LogMessageType type) const
{
  return m_logColors[type];
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void LoggingList::setMaxLogLines(unsigned int maxLogLines)
{
  m_maxLogLines = maxLogLines;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned int LoggingList::getMaxLogLines() const
{
  return m_maxLogLines;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

unsigned int LoggingList::getLogLinesCounter() const
{
  return m_logLinesCounter;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void LoggingList::clearLog()
{
  m_logLinesCounter = 0;
  this->clear();
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void LoggingList::appendToLog(const QString & string, LogMessageType type)
{
  // get date and time
  QString date = QDate::currentDate().toString("MM/dd/yyyy");
  QString time = QTime::currentTime().toString("hh:mm:ss");
  
  // copy useful value to avoid repeated accesses
  // to the same elements in the hash maps
  QString color = m_logColors[type].name();
  QString prefix = m_prefixes[type];
  
  // split lines to a string list (all lines must be begin with 
  // date, time and prefix)
  QStringList list = string.split("\n", QString::SkipEmptyParts);
  
  for(int i = 0 ; i < list.size() ; i++)
  {
    QString str;
    QString text;
    
    if(m_maxLogLines != 0)
    {
      // if log has 100000 lines, we clear it (to save memory)
      if(m_logLinesCounter == m_maxLogLines)
        this->clearLog();
      
      else
        m_logLinesCounter++;
    }
    
    // build the string to append
    text = prefix + " " + list.at(i);
    
    str = QString("<font color=\"%1\">%2</font>").arg(color).arg(text);
    str.prepend(QString("[%1 %2] ").arg(date).arg(time));
    
    this->append(QString("<font face=\"monospace\">%3</font>").arg(str));
  }
}

/****************************************************************************
 
 PUBLIC SLOTS
 
 ****************************************************************************/

void LoggingList::message(const QString & message, bool fromServer)
{
  LogMessageType type = TYPE_NORMAL;
  QString prefix;
  QString messageToAppend = message;
  
  if(fromServer)
  {
    if(message.startsWith("NOTICE"))
    {
      type = TYPE_SERVER_NOTICE;
      messageToAppend.remove(QRegExp("^NOTICE"));
    }
    
    else if(message.startsWith("INFO"))
    {
      type = TYPE_SERVER_INFO;
      messageToAppend.remove(QRegExp("^INFO"));
    }
    
    else if(message.startsWith("WARN"))
    {
      type = TYPE_SERVER_WARN;
      messageToAppend.remove(QRegExp("^WARN"));
    }
    
    else
      type = TYPE_SERVER_NOTICE;
  }
  
  this->appendToLog(messageToAppend, type);
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void LoggingList::error(const QString & message, bool fromServer)
{
  LogMessageType type = fromServer ? TYPE_SERVER_ERROR : TYPE_ERROR;
  
  this->appendToLog(message, type);
}
