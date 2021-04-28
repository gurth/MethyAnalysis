#include "mainwindow.h"
#include "my_config.h"
#include "interface.h"
#include "../tool/MethyProfile/bed.h"
#include <QMessageBox>
#include <cstring>

using namespace bed;

 m_Interface::m_Interface()
 {

 }

 m_Interface::m_Interface(const char* bedname, const char* gff3name, const char* outputname, bool have_promoter
                          , size_t length_promoter, bool do_single_analysis, const char* sglist)
 {
    set_arg(bedname,gff3name,outputname,have_promoter,length_promoter,do_single_analysis,sglist);
 }

 void m_Interface::set_arg(const char* bedname, const char* gff3name, const char* outputname, bool have_promoter
                           , size_t length_promoter, bool do_single_analysis, const char* sglist)
 {
     strcpy(this->bedname, bedname);
     strcpy(this->gff3name, gff3name);
     strcpy(this->outputname, outputname);
     this->have_promoter=have_promoter;
     this->length_promoter=length_promoter;
     this->do_single_analysis=do_single_analysis;
     strcpy(this->sglist, sglist);

     pThis=this;
 }

 m_Interface::~m_Interface()
 {

 }

void m_Interface::run()
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

    bedfile.process(nullptr, Method::tag);

    sprintf(msg_buff, "Tagging finished.\n\nTotal cost: %.4lf ms\n", bedfile.latest_time_cost);

    emit send_msgBox(msg_buff);

    bedfile.savechrList();

    bedfile.process(gff3name, outputname, Method::profile);

    sprintf(msg_buff, "Profile generation finished.\n\nTotal cost: %.4lf ms\n", bedfile.latest_time_cost);

    emit send_msgBox(msg_buff);
}

void m_Interface::m_setProgress(double progress)
{
    emit pThis->send_setProgress(progress);
}

void m_Interface::m_initProgress(double progress, const char *info)
{
    emit pThis->send_setProgress(progress);
    emit pThis->send_setMsg(info);
}

void m_Interface::m_errorExit(int m_error)
{
    emit pThis->send_error(m_error);
}
