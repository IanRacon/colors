#include "Configuration.h"

const std::string section = "colors";

Configuration& Configuration::getInstance()
{
    static Configuration instance("configuration.ini");
    return instance;
}

Configuration::Configuration(const std::string& filename):reader(filename)
{
}
int Configuration::getInt(const std::string& name) const 
{
    return reader.GetInteger(section, name, 0);
}
double Configuration::getDouble(const std::string& name) const 
{
    return reader.GetReal(section, name, 0.0);
}
std::string Configuration::getString(const std::string& name) const 
{
    return reader.Get(section, name, "");
}
bool Configuration::getBool(const std::string& name) const 
{
    return reader.GetBoolean(section, name, true);
}