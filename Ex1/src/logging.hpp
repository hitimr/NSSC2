#include <iostream>
#include <string>
#include <fstream>
#include <time.h>

class Logger
{
public:
    int m_rank;
    bool master;
    std::string m_separator = ";";
    std::string m_logFile;
    std::ofstream m_outStream;

    Logger(int rank);
    ~Logger();

    const std::string currentDateTime();

    void add(
        std::string res, 
        std::string nprocs, 
        std::string runtimeval, 
        std::string errorval,
        std::string errorMax,
        std::string residualNorm,
        std::string residualMax,
        std::string average_iteration_time,
        std::string iterations,
        std::string dtype,
        std::string mpi)
    {
        if(master)
        {
            m_outStream.open(m_logFile, fstream::app);
            m_outStream 
                << res << ";" 
                << nprocs << ";" 
                << runtimeval << ";" 
                << errorval << ";" 
                << errorMax << ";" 
                << residualNorm << ";" 
                << residualMax << ";" 
                << average_iteration_time << ";" 
                << iterations << ";" 
                << dtype << ";" 
                << mpi
                << std::endl;
            m_outStream.close();
        }
};

Logger::Logger(int rank)
{
    m_rank = rank;
    if(m_rank == 0) master = true;
    else master = false;

    //m_logFile = string("./out/logs/log_" + currentDateTime() + ".csv");
    m_logFile = "./out/logs/log.csv";

    if(master)
    {
        m_outStream.open(m_logFile, fstream::app);
    }
}

Logger::~Logger()
{
    m_outStream.close();
}

// Get current date/time, format is YYYY-MM-DD.HH:mm:ss
const std::string Logger::currentDateTime() {
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    strftime(buf, sizeof(buf), "%Y-%m-%d-%H-%M-%S", &tstruct);

    return buf;
}

