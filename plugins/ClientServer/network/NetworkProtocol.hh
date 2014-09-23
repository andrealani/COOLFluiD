#ifndef COOLFluiD_network_NetworkProtocol_h
#define COOLFluiD_network_NetworkProtocol_h

//////////////////////////////////////////////////////////////////////////////

#include "Config/BuilderParserRules.hh"

#include "ClientServer/network/NetworkAPI.hh"
#include "ClientServer/network/NetworkFrameType.hh"

class QString;

namespace COOLFluiD
{
  namespace network
  {
    
    //////////////////////////////////////////////////////////////////////////////
    
    /// @brief Defines the network m_protocol rules
    
    /// @author Quentin Gasper
    
    class Network_API NetworkProtocol : public COOLFluiD::Config::BuilderParserRules
    {
    public:
      
      /// @brief Constructor
      NetworkProtocol();
      
      /// @brief Convert an unsigned int to the corresponding type.
      
      /// @param type The integer to convert
      /// @return Returns the converted type, or @c #NETWORK_NO_TYPE if the type
      /// could not be converted.
      NetworkFrameType convertToType(unsigned int m_type) const;
      
      /// @brief Convert a string to the corresponding type.
      
      /// @param type The string to convert
      /// @return Returns the converted type, or @c #NETWORK_NO_TYPE if the type
      /// could not be converted.
      NetworkFrameType convertToType(const QString & typeName) const;
      
    private:
      
      /// @brief Builds the rules
      void buildRules();
      
    };
    
    //////////////////////////////////////////////////////////////////////////////
    
  } // namespace network
} // namespace COOLFluiD

//////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_network_NetworkProtocol_h