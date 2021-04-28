#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QProgressBar>
#include <QLabel>

#include "interface.h"

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
    /*
    static void MethyProfileInterface(const char* bedname, const char* gff3name, const char* outputname, bool have_promoter
                               , size_t length_promoter, bool do_single_analysis, const char* sglist);

    static int infoResult(const char* info);*/

private slots:
    void rev_setProgress(int progress);

    void rev_setMsg(const char* info);

    void rev_msgBox(const char* info);

    void rev_error(int m_error);

private:
    Ui::MainWindow *ui;
    QWidget *settingWidget = nullptr;
    QProgressBar *m_progressbar = nullptr;
    QLabel *m_labelstatus = nullptr;

    m_Interface mInf;

};

#endif // MAINWINDOW_H
