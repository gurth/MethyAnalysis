#ifndef INTERFACE_H
#define INTERFACE_H

#include <QThread>
#include <limits.h>

#ifdef _UNIX_PLATFORM_

#define BED_MAX_PATH NAME_MAX

#elif defined(_WIN32_PLATFORM_)

# include <Windows.h>
#define BED_MAX_PATH MAX_PATH

#endif //!_UNIX_PLATFORM_

class m_Interface:public QThread
{
    Q_OBJECT
private:
    char bedname[BED_MAX_PATH]={0};
    char gff3name[BED_MAX_PATH]={0};
    char outputname[BED_MAX_PATH]={0};
    bool have_promoter=false;
    size_t length_promoter=0;
    bool do_single_analysis=false;
    char sglist[BED_MAX_PATH]={0};

    static m_Interface* pThis;
public:
    m_Interface();
    m_Interface(const char* bedname, const char* gff3name, const char* outputname, bool have_promoter
                , size_t length_promoter, bool do_single_analysis, const char* sglist);
    ~m_Interface();
    void run() override;

    void set_arg(const char* bedname, const char* gff3name, const char* outputname, bool have_promoter
                  , size_t length_promoter, bool do_single_analysis, const char* sglist);
    static void m_setProgress(double progress);
    static void m_initProgress(double progress, const char* info);
    static void m_errorExit(int m_error);
signals:
    void send_setProgress(int progress);
    void send_setMsg(const char* info);
    void send_msgBox(const char* info);
    void send_error(int m_error);
};

#endif // INTERFACE_H
