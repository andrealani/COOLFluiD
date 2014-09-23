#ifndef COOLFluiD_server_Worker_h
#define COOLFluiD_server_Worker_h

#include <QHash>
#include "server/WorkerThre.hh"

class QString;

////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD
{
  namespace server
  {
    
    class WorkerThread;
    
    ////////////////////////////////////////////////////////////////////////////
    
    class Worker 
    {    
    public:
      
      Worker();
      
      QString getString(const QString & string = QString());
      
      void simulate(MgrWorkerProtocol & rules);
      
      /// @brief Gives the current value.
      
      /// @return Returns the current value.
      int getNumber() ;
      
      /// @brief Gives the sum of all thread numbers.
      
      /// @return Returns the sum of all thread numbers.
      int getSum();
      
      QString getHeader() const;
      
      void openPort(const QString & remoteRole);
      
      void connectToPort(const QString & portName, const QString & remoteRole);
      
      void setRole(const QString & header);
      
      QString getRole() const;
      
    private:
      
      /// @brief Current value
      int m_number;
      
      /// @brief Sum of all thread values
      int m_sum;
      
      /// @brief Process rank
      int m_rank;
      
      MPI::Intercomm m_parent;
      
      QHash<QString, MPI::Intercomm> m_global;
      
      int m_count;
      
      char m_string[32768]; 
      
      QString m_role;
      
      
    };
    
    ////////////////////////////////////////////////////////////////////////////
    
  } // namespace server
} // namespace COOLFluiD

////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_server_Worker_h
