#pragma once
#include <string>
#include "INIReader.h"
class Configuration{
public:
    static Configuration& getInstance();
    int getInt(const std::string& name) const;
    bool getBool(const std::string& name) const;
    double getDouble(const std::string& name) const;
    std::string getString(const std::string& name) const;
private:
    Configuration(const std::string& fileName);
    INIReader reader;
};