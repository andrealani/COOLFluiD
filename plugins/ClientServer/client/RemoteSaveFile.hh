#ifndef COOLFluiD_client_RemoteSaveFile_h
#define COOLFluiD_client_RemoteSaveFile_h

#include "ClientServer/client/RemoteFSBrowser.hh"

class QModelIndex;

/////////////////////////////////////////////////////////////////////////////

namespace COOLFluiD
{
  namespace client
  {
    class TypeAndNameDialog;
    class ClientNetworkComm;
    /////////////////////////////////////////////////////////////////////////////
    
    /// @brief This class is a dialog that allows user to to select a remote
    /// place to save a file.
    
    /// This class subclasses @c RemoteFSBrowser. User can create new directories,
    /// enter a name for the new file or choose to overwrite a file. This class
    /// does not allow user to call @c #setIncludeFiles, @c #setIncludeNoExtension
    /// and @c #showMultipleSelect base class methods. 
    
    /// @author Quentin Gasper
    class RemoteSaveFile : public RemoteFSBrowser
    {
      Q_OBJECT
      
    public:
      
      /// @brief Constructor
      
      /// @param parent Parent window. May be @c NULL.
      RemoteSaveFile(const QModelIndex & index, QMainWindow * parent = NULL);
      
      /// @brief Destructor
      
      /// Frees all allocated memory. Parent is not destroyed.
      ~RemoteSaveFile();
      
    protected:
      
      /// @brief Checks if the selection is valid. 
      
      /// This method overrides base class method. Three cases are possible:
      /// @li if the given name is a directory and user has entered a file name,
      /// the file name is appended to the directory path and the method returns
      /// @c #POLICY_VALID.
      /// @li if the given name is a directory but user has not entered a file 
      /// name, the method returns @c #POLICY_ENTER_DIRECTORY.
      /// @li if the given name is a file: a confirmation to overwrite the is asked.
      /// If user confirms, @c #POLICY_VALID is returned; otherwise 
      /// @c #POLICY_NOT_VALID is returned.
      /// @note If the given name is a file, any file name entered by user is 
      /// ignored.
      /// @param name Name to check
      /// @param isDir If @c true, @c name is a directory; otherwise, it is a file.
      /// @return Returns the validation as described above.
      virtual ValidationPolicy isAcceptable(const QString & name, bool isDir);
      
      /// @brief Reinitializes internal data to their default value.
      
      /// This method overrides base class method.
      virtual void reinitValues();
      
      /// @brief Give the selected file.
      
      /// This method overrides base class method.
      /// @return Returns the selected file, or an empty string if the last call
      /// to @c #isAcceptable did not return @c #POLICY_VALID or if this method was
      /// never called.
      virtual QString getSelectedFile() const;
      
      private slots:
      
      /// @brief Asks user to enter a file name and select file type.
      
      /// File type associated extention is appended to the name if needed. If
      /// user enters a name that already exists, he has to confirm his entry.
      void btFileNameClick();
      
      /// @brief Asks user to enter a new directory name.
      
      /// If a name is entered, a request to create that directory is sent to 
      /// the server.
      void btNewDirectoryClicked();
      
    private:
      
      /// @brief "Set file name" button
      QPushButton * m_btFileName;
      
      /// @brief "New directory" button
      QPushButton * m_btNewDirectory;
      
      /// @brief File name entered by user.
      QString m_fileName;
      
      /// @brief Selected file and its path.
      QString m_selectedFile;
      
      /// @brief Dialog used to ask the new file name
      TypeAndNameDialog * m_fileNameDialog;
      
      
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
  }
}  

/////////////////////////////////////////////////////////////////////////////


#endif // COOLFluiD_client_RemoteSaveFile_h
