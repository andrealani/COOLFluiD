#ifndef COOLFluid_client_FilesPanel_h
#define COOLFluid_client_FilesPanel_h

/////////////////////////////////////////////////////////////////////////////

#include <QWidget>
#include <QHash>

class QComboBox;
class QGridLayout;
class QHBoxLayout;
class QLabel;
class QListView;
class QPushButton;
class QStringList;
class QStringListModel;
class QVBoxLayout;

namespace COOLFluiD
{
  namespace client
  {
    
    class RemoteOpenFile;
    
    enum ComboItemType
    {
      ITEM_ADD_ITEMS,
      ITEM_REMOVE_ITEMS,
      ITEM_INVERT_SELECTION,
      ITEM_CLEAR_SELECTION,
      ITEM_CLEAR_LIST
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
    /// @brief This class represents a graphical component that allows user to 
    /// maintain a list of files.
    
    /// This component is composed of a list m_view (to display the list) and four
    /// m_buttons : one to add files to the list, one to remove files from the list,
    /// one to clear the list and the last one to clear list selection (selected
    /// m_items become unselected).@n @n
    /// To add files, the component uses @c RemoteFSBrowser class and allows the 
    /// dialog to select several files at the same time. The list m_view also 
    /// allows the user to select several files at the same time.
    /// This class subclasses @c QWidget so it can be used as any other graphical
    /// component.
    
    class FilesPanel : public QWidget
    {
      Q_OBJECT
      
    public:
      /// @brief Constructor
      
      /// @param includeFiles If @c true, files are requested.
      /// @param extensions List of wanted files extension. May be empty if all 
      /// files with extension are wanted.
      /// @param includeNoExtension If @c true, files without extension are 
      /// requested too. The parameter has no effect if @c includeFiles is 
      /// @c false.
      /// @param parent Parent widget. May be @c NULL.
      FilesPanel(bool includeFiles, const QStringList & extensions, 
                 bool includeNoExtension, QWidget * parent = NULL);
      
      /// @brief Destructor
      
      /// Frees all alocated memory. The parent is not destroyed
      ~FilesPanel();
      
      /// @brief Gives the list of files.
      /// @return Returns the list of files.
      virtual QStringList getFilesList() const;
      
      /// @brief Changes the files list
      
      /// The existing files list is deleted.
      /// @param filesList The new files list.
      void setFilesList(const QStringList & filesList);
      
    protected:
      /// @brief Allows a subclass to change m_buttons text suffix.
      
      /// "Add" and "Remove" m_buttons are affected. For example, if @c name is
      /// "files", these button text will become respectively "Add files" and 
      /// "Remove files". The m_buttons behaviour remains the same (add/remove m_items
      /// from the list).
      /// @param name Text suffix. May be empty.
      void setButtonNames(const QString & name);
      
      /// @brief Add files to the list.
      void addFile();
      
      /// @brief Removes selected item from the list.
      void removeFile();
      
      /// @brief Clears the list.
      void clearList();
      
      /// @brief Inverts the current selection.
      void invertSelection();
      
      private slots:
      
      /// @brief Slot called when user clicks on "Ok" button.
      void btOkClicked();
      
    private:
      
      /// @brief The dialog to select the files // 
      RemoteOpenFile * openFileDialog;
      
      /// @brief List m_view that displays the files 
      QListView * m_filesListView;
      
      /// @brief The model that maintains the files list.
      
      /// This is the model of the m_view. 
      QStringListModel * m_filesListModel;
      
      /// @brief Main m_layout 
      QGridLayout * m_mainLayout;
      
      /// @brief Layout that contains the comboBox and "Ok" button 
      QHBoxLayout * m_buttonsLayout;
      
      QVBoxLayout * m_actionsLayout;
      
      /// @brief Widget that displays @c m_buttonsLayout 
      QWidget * m_buttonsWidget;
      
      /// @brief ComboBox with actions 
      QComboBox * m_comboActions;
      
      /// @brief "OK" butons 
      QPushButton * m_btOk;
      
      /// @brief "Actions" label 
      QLabel * m_labActions;
      
      /// @brief ComboBox m_items 
      QHash<int, QString> m_comboActionItems;
      
      
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
  } // namespace client
} // namespace COOLFluiD

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluid_client_FilesPanel_h
