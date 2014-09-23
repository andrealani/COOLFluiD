#ifndef COOLFluiD_network_HostInfos_h
#define COOLFluiD_network_HostInfos_h

//////////////////////////////////////////////////////////////////////////////

class QString;

namespace COOLFluiD
{
  namespace network
  {
    
    //////////////////////////////////////////////////////////////////////////////
    
    /// @brief Holds host information.
    
    /// @author Quentin Gasper
    
    struct HostInfos
    {
    public:
      
      /// @brief Hostname
      QString m_hostname;
      
      /// @brief Number of slots the host has.
      unsigned int m_nbSlots;
      
      /// @brief Maximum number of slot that can be allocated
      unsigned int m_maxSlots;
      
      /// @brief Constructor
      HostInfos(const QString & hostname = QString(), int nbSlots = 0, 
                int maxSlots = 0);
      
      //////////////////////////////////////////////////////////////////////////////
      
    }; // struct HostInfos
    
    //////////////////////////////////////////////////////////////////////////////
    
  } // namespace network
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_network_HostInfos_h