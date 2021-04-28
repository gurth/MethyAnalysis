#include "mainwindow.h"
#include "./ui_mainwindow.h"
#include <QFileDialog>
#include <QMessageBox>
#include "my_config.h"
#include "my_error.hpp"

m_Interface* m_Interface::pThis=nullptr;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    settingWidget = ui->m_settingWidget;

    init();
}

void MainWindow::init()
{
    settingWidget->hide();
    this->resize(this->width(),0);

    m_progressbar=ui->m_progressBar;
    m_labelstatus=ui->m_labelStatus;

    connect(&mInf, &m_Interface::send_error, this, &MainWindow::rev_error);

    connect(&mInf, &m_Interface::send_setProgress, this, &MainWindow::rev_setProgress);
    connect(&mInf, &m_Interface::send_setMsg, this, &MainWindow::rev_setMsg);
    connect(&mInf, &m_Interface::send_msgBox, this, &MainWindow::rev_msgBox);
}

MainWindow::~MainWindow()
{
    mInf.wait();
    delete ui;
}

void MainWindow::on_m_buttonShowSetting_clicked()
{
    if(settingWidget->isHidden())
    {
        settingWidget->show();
        ui->m_widgetTitle->hide();
        ui->m_buttonShowSetting->setText("Setting <<<<");
    }
    else
    {
        settingWidget->hide();
        ui->m_widgetTitle->show();
        ui->m_buttonShowSetting->setText("Setting >>>>");
    }
    this->resize(this->width(),0);
}

void MainWindow::on_m_pushButtonBED_clicked()
{
    QString bedfileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                        "",
                                                        tr("BED File (*.bed *.txt) ;; (*.*)"));
    if(!bedfileName.isEmpty())
    {
        if(ui->m_comboBoxBED->findText(bedfileName) == -1)
            ui->m_comboBoxBED->addItem(bedfileName);
        ui->m_comboBoxBED->setCurrentText(bedfileName);
        if(ui->m_comboBoxOUT->currentIndex()==-1)
        {
            bedfileName+=".methyprofile.txt";
           if(ui->m_comboBoxOUT->findText(bedfileName) == -1)
               ui->m_comboBoxOUT->addItem(bedfileName);
           ui->m_comboBoxOUT->setCurrentText(bedfileName);
        }
    }
}

void MainWindow::on_m_pushButtonGFF3_clicked()
{
    QString gff3fileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                        "",
                                                        tr("GFF3 File (*.gff3 *.txt) ;; (*.*)"));
    if(!gff3fileName.isEmpty())
    {
        if(ui->m_comboBoxGFF3->findText(gff3fileName) == -1)
            ui->m_comboBoxGFF3->addItem(gff3fileName);
        ui->m_comboBoxGFF3->setCurrentText(gff3fileName);
    }
}

void MainWindow::on_m_pushButtonOUT_clicked()
{
    QString outfileName = QFileDialog::getSaveFileName(this, tr("Open File"),
                                                       "",
                                                       tr("MethyProfile (*.txt)"));
    if(!outfileName.isEmpty())
    {
        if(ui->m_comboBoxOUT->findText(outfileName) == -1)
            ui->m_comboBoxOUT->addItem(outfileName);
        ui->m_comboBoxOUT->setCurrentText(outfileName);
    }
}


void MainWindow::on_m_pushButtonSGFile_clicked()
{
    QString listfileName = QFileDialog::getOpenFileName(this, tr("Open File"),
                                                        "",
                                                        tr("Gene list File (*.list *.txt) ;; (*.*)"));
    if(!listfileName.isEmpty())
    {
        if(ui->m_comboBoxSGFile->findText(listfileName) == -1)
            ui->m_comboBoxSGFile->addItem(listfileName);
        ui->m_comboBoxSGFile->setCurrentText(listfileName);
    }
}

void MainWindow::on_m_pushButtonStart_clicked()
{
    QString bedfileName =  ui->m_comboBoxBED->currentText();
    QString gff3fileName =  ui->m_comboBoxGFF3->currentText();
    QString outfileName =  ui->m_comboBoxOUT->currentText();

    if(bedfileName.isEmpty())
    {
        QMessageBox::critical(NULL,  "Error",  "BED file name cannot be empty.", QMessageBox::Yes, QMessageBox::Yes);
        return;
    }

    if(gff3fileName.isEmpty())
    {
        QMessageBox::critical(NULL,  "Error",  "GFF3 file name cannot be empty.", QMessageBox::Yes, QMessageBox::Yes);
        return;
    }

    bool have_promoter=ui->m_checkBoxPromoter->isChecked();
    bool do_single_analysis=ui->m_checkBoxSG->isChecked();
    size_t length_promoter=0;
    QString listfileName = "";

    if(have_promoter)
    {
        length_promoter=ui->m_lineEditPromoterLen->text().toULong();
    }
    if(do_single_analysis)
    {
        listfileName=ui->m_comboBoxSGFile->currentText();
        if(listfileName.isEmpty())
        {
            QMessageBox::critical(NULL,  "Error",  "Single gene list file name cannot be empty.", QMessageBox::Yes, QMessageBox::Yes);
            return;
        }
    }

    try
    {
        mInf.set_arg(
                    bedfileName.toLocal8Bit().data(),
                    gff3fileName.toLocal8Bit().data(),
                    outfileName.toLocal8Bit().data(),
                    have_promoter,
                    length_promoter,
                    do_single_analysis,
                    listfileName.toLocal8Bit().data());
        mInf.start();
    }
    catch (const e_ErrorCode& ex)
    {
        QMessageBox::critical(NULL,  "Error",  ex.what(), QMessageBox::Yes, QMessageBox::Yes);
        return;
    }

}

void MainWindow::rev_setMsg(const char *info)
{
    ui->m_labelStatus->setText(info);
}

void MainWindow::rev_msgBox(const char *info)
{
    QMessageBox message(QMessageBox::NoIcon,  "Result",  info);
    message.setIconPixmap(QPixmap(":/res/methyprofile_64.png"));
    int rsl=message.exec();
}

void MainWindow::rev_error(int m_error)
{
    char msg_buff[MAX_INFO]={0};
    sprintf(msg_buff, "Unexpect exit with error code %d. \n", m_error);
    QMessageBox::critical(NULL,  "Error",  msg_buff, QMessageBox::Yes, QMessageBox::Yes);
    switch (m_error)
    {
        default:
            break;
    }
    throw e_ErrorCode(m_error);
}

void MainWindow::rev_setProgress(int progress)
{
    ui->m_progressBar->setValue(progress);
}
