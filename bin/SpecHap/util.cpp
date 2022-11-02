#include "util.h"
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

void logging(std::ostream &stream, const std::string & message)
{
    std::time_t t = std::time(nullptr);

    stream << "[SpecHap " << std::put_time(std::localtime(&t), "%Y:%m:%d %H:%M:%S]") << message << std::endl;
}
