#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QProgressBar>
#include <QLabel>

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void on_m_pushButtonBED_clicked();

    void on_m_buttonShowSetting_clicked();

    void on_m_pushButtonGFF3_clicked();

    void on_m_pushButtonOUT_clicked();

    void on_m_pushButtonSGFile_clicked();

    void on_m_pushButtonStart_clicked();

private:
    void init();
    void MethyProfileInterface(const char* bedname, const char* gff3name, const char* outputname, bool have_promoter
                               , size_t length_promoter, bool do_single_analysis, const char* sglist);
    static void m_setProgress(double progress);
    static void m_initProgress(double progress, const char* info);
    static void m_errorExit(int m_error);

    static int infoResult(const char* info);

private:
    Ui::MainWindow *ui;
    QWidget *settingWidget = nullptr;
    QProgressBar *m_progressbar = nullptr;
    QLabel *m_labelstatus = nullptr;

    static MainWindow* pThis;
};

#endif // MAINWINDOW_H
