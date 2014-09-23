#ifndef COOLFluiD_client_CloseConfirmationPanel_h
#define COOLFluiD_client_CloseConfirmationPanel_h

/////////////////////////////////////////////////////////////////////////////

#include <QWidget>
#include <QVector>
#include <QDialogButtonBox>

class QLabel;
class QPushButton;
class QVBoxLayout;
class QGridLayout;

namespace COOLFluiD
{
  namespace client
  {
    
    struct CloseConfirmationInfos;
    
    /////////////////////////////////////////////////////////////////////////////
    
    /// @brief This class is a base class for confirmation m_panels.
    
    /// This class contains pure virtual methods and can not be instancied
    /// directly. Please, use @c CommitConfirmationPanel, 
    /// @c SaveConfirmationPanel and @c ShutdownConfirmationPanel subclasses
    /// instead.
    class CloseConfirmationPanel : public QWidget
    {
      Q_OBJECT
      
    public:
      
      /// @brief Constructor
      
      /// @param title Panel title
      /// @param parent Parent dialog. May be @c NULL.
      CloseConfirmationPanel(const QString title, QDialog * parent = 0);
      
      /// @brief Destructor
      
      /// Frees all allocated memory. Parent is not destroyed.
      ~CloseConfirmationPanel();
      
      /// @brief Indicates whether user correctly answered the question.
      
      /// @return Returns @c true if user correctly answered the question; 
      /// otherwise returns @c false.
      virtual bool isAccepted() const = 0;
      
      /// 
      virtual void getData(CloseConfirmationInfos & infos) const = 0;
      
      virtual void setData(const CloseConfirmationInfos & infos) = 0;
      
    protected:
      
      void setText(const QString & text);
      
      void setHelp(const QString & helpText);
      
      void addButton(const QString & text, QDialogButtonBox::ButtonRole role,
                     const char * slot);
      
      void addWidget(QWidget * widget);
      
      void setWidgetBesideButtons(QWidget * widget);
      
      // does not hide widgetBesideButtons
      virtual void hideComponents(bool hide);
      
      private slots:
      
      void help();
      
      void titleClicked(const QString & link);
      
    signals:
      
      void showHelp(const QString & m_helpText);
      
      void resized();
      
    private:
      
      /// @brief Clickable label for the panel title
      QLabel * m_labTitle;
      
      /// @brief Label for some text that explains the problem
      QLabel * m_labText;
      
      /// @brief Help text
      QString m_helpText;
      
      /// @brief Main m_layout
      QVBoxLayout * m_layout;
      
      /// @brief Layout for m_buttons
      QGridLayout * m_buttonsLayout;
      
      /// @brief Panel m_buttons
      QDialogButtonBox * m_buttons;
      
      /// @brief Widget that appears beside (on the left) panel m_buttons.
      QWidget * m_besideButtonsWidget;
      
      /// @brief Widget for m_buttons
      
      /// Its m_layout is @c m_buttonsLayout
      QWidget * m_buttonsWidget;
      
      /// @brief Vector of m_buttons that were added by the subclasse
      QVector<QPushButton *> m_buttonsVector;
      
      /// @brief Indicates whether panel components are m_hidden.
      bool m_hidden;
      
    };
    
    /////////////////////////////////////////////////////////////////////////////
    
  }
}

/////////////////////////////////////////////////////////////////////////////

#endif // COOLFluiD_client_CloseConfirmationPanel_h
