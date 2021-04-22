#include "mainwindow.h"
#include "my_config.h"
#include "../tool/MethyProfile/bed.h"
#include <QMessageBox>
#include <cstring>

using namespace bed;

void MainWindow::m_setProgress(double progress)
{
    pThis->m_progressbar->setValue(progress);
    pThis->m_progressbar->update();
}

void MainWindow::m_initProgress(double progress, const char *info)
{
    pThis->m_progressbar->setValue(progress);
    pThis->m_progressbar->update();
    QString str_status=info;
    pThis->m_labelstatus->setText(str_status);
    pThis->m_labelstatus->update();
}

void MainWindow::m_errorExit(int m_error)
{
    char msg_buff[MAX_INFO]={0};
    sprintf(msg_buff, "Unexpect exit with error code %d. \n", m_error);
    QMessageBox::critical(NULL,  "Error",  msg_buff, QMessageBox::Yes, QMessageBox::Yes);
    switch (m_error)
    {
        default:
            break;
    }
    exit(m_error);
}


int MainWindow::infoResult(const char *info)
{
    QMessageBox message(QMessageBox::NoIcon,  "Result",  info);
    message.setIconPixmap(QPixmap(":/res/methyprofile_64.png"));
    return message.exec();
}

void MainWindow::MethyProfileInterface(const char* bedname, const char* gff3name, const char* outputname, bool have_promoter
                           , size_t length_promoter, bool do_single_analysis, const char* sglist)
{
    BED bedfile(m_errorExit);
    char msg_buff[MAX_INFO]={0};

    bedfile.have_promoter=have_promoter;
    if(have_promoter)
    {
        bedfile.promoterLen=length_promoter;
    }
    bedfile.do_single_analyse=do_single_analysis;
    if(do_single_analysis)
    {
        bedfile.loadSingleList(sglist);
    }
    bedfile.bedfileOpen(bedname);

    bedfile.initProgress=m_initProgress;
    bedfile.setProgress=m_setProgress;

    pThis = this;
    bedfile.process(nullptr, Method::tag);

    sprintf(msg_buff, "Tagging finished.\n\nTotal cost: %.4lf ms\n", bedfile.latest_time_cost);

    infoResult(msg_buff);

    bedfile.savechrList();
    bedfile.process(gff3name, outputname, Method::profile);

    sprintf(msg_buff, "Profile generation finished.\n\nTotal cost: %.4lf ms\n", bedfile.latest_time_cost);

    infoResult(msg_buff);

    return;
}
