#ifndef __LOG1_H__
#define __LOG1_H__

#include <sstream>
#include <string>
#include <stdio.h>
#include <iostream>
using std::cout;
using std::endl;

inline std::string NowTime();

enum TLogLevel {logNONE=0,logERROR=1, logWARN=2, logINFO=3, logDEBUG=4, logDEBUG1=5, logDEBUG2=6, logDEBUG3=7, logDEBUG4=8};

class Log
{
public:
    Log();
    virtual ~Log();
    std::ostringstream& Get(TLogLevel level = logINFO);
public:
    static TLogLevel& ReportingLevel();
    static std::string ToString(TLogLevel level);
    static TLogLevel FromString(const std::string& level);
    static TLogLevel FromInt(const int level);
protected:
    std::ostringstream os;
private:
    Log(const Log&);
    Log& operator =(const Log&);
};

inline Log::Log()
{
}

inline std::ostringstream& Log::Get(TLogLevel level)
{
  os <<"- " << NowTime();
  os << " " << ToString(level) << ": ";
  os << std::string(level > logDEBUG ? level - logDEBUG : 0, ' ');
  return os;
}

inline Log::~Log()
{
  //os << std::endl;
  cout<<os.str();//<<endl;
}

inline TLogLevel& Log::ReportingLevel()
{
    static TLogLevel reportingLevel = logDEBUG4;
    return reportingLevel;
}

inline std::string Log::ToString(TLogLevel level)
{
  static const char* const buffer[] = {"NONE","ERROR", "WARN", "INFO", "DEBUG", "DEBUG1", "DEBUG2", "DEBUG3", "DEBUG4"};
  return buffer[level];
}

inline TLogLevel Log::FromString(const std::string& level)
{
    if (level == "DEBUG4")
        return logDEBUG4;
    if (level == "DEBUG3")
        return logDEBUG3;
    if (level == "DEBUG2")
        return logDEBUG2;
    if (level == "DEBUG1")
        return logDEBUG1;
    if (level == "DEBUG")
        return logDEBUG;
    if (level == "INFO")
        return logINFO;
    if (level == "WARN")
        return logWARN;
    if (level == "ERROR")
        return logERROR;
    if (level == "NONE")
        return logNONE;
    Log().Get(logWARN) << "Unknown logging level '" << level << "'. Using INFO level as default.";
    return logINFO;
}

inline TLogLevel Log::FromInt(const int level)
{
  if (level == 8)
    return logDEBUG4;
  if (level == 7)
    return logDEBUG3;
  if (level == 6)
    return logDEBUG2;
  if (level == 5)
    return logDEBUG1;
  if (level == 4)
    return logDEBUG;
  if (level == 3)
    return logINFO;
  if (level == 2)
    return logWARN;
  if (level == 1)
    return logERROR;
  if (level == 0)
    return logNONE;
  Log().Get(logWARN) << "Unknown logging level '" << level << "'. Using INFO level as default.";
  return logINFO;
}

typedef Log FILELog;

#define FILE_LOG(level) \
  if (level > FILELog::ReportingLevel()) ;       \
  else Log().Get(level)

#include <sys/time.h>

inline std::string NowTime()
{
    char buffer[11];
    time_t t;
    time(&t);
    tm r = {0};
    strftime(buffer, sizeof(buffer), "%X", localtime_r(&t, &r));
    struct timeval tv;
    gettimeofday(&tv, 0);
    char result[100] = {0};
    std::sprintf(result, "%s.%03ld", buffer, (long)tv.tv_usec / 1000); 
    return result;
}

#endif //__LOG_H__
