#ifndef METHYANALYSIS_MY_ERROR_H
#define METHYANALYSIS_MY_ERROR_H

#include<exception>
#include<string>

class e_ErrorCode : public std::exception
{
private:
    std::string m_msg;
public:
    explicit e_ErrorCode(int m_error)
    {
        this->m_msg="Process exit with code ";
        this->m_msg+=std::to_string(m_error);
        this->m_msg+=". ";
    }
    const char* what() const noexcept override
    {
        return m_msg.c_str();
    }
};

#endif //!METHYANALYSIS_MY_ERROR_H
