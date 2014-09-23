#ifndef COOLFluid_client_LibrariesPanel_h
#define COOLFluid_client_LibrariesPanel_h

/////////////////////////////////////////////////////////////////////////////

#include "ClientServer/client/FilesPanel.hh"

class QStringList;

namespace COOLFluiD
{
  namespace client
  {
    
    /////////////////////////////////////////////////////////////////////////////
    
    /// @brief This class represents a graphical component that allows user to 
    /// maintain a list of libraries.
    
    class LibrariesPanel : public FilesPanel
    {
    public:
      /// @brief Constructor
      
      /// @param parent Parent widget. May be @c NULL.
      LibrariesPanel(QWidget * parent = NULL);
      
      /// @brief Destructor
      
      /// Frees all alocated memory. The parent is not destroyed
      ~LibrariesPanel();
      
      /// @brief Gives the list of files.
      /// @return Returns the list of files.
      virtual QStringList getFilesList() const;
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
  } // namespace client
} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluid_client_LibrariesPanel_h
