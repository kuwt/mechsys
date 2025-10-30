#include "json.hpp"
#include <iostream>

template < class T >
std::ostream& operator << (std::ostream& os, const std::vector<T>& v) 
{
    os << "[";
    for (typename std::vector<T>::const_iterator ii = v.begin(); ii != v.end(); ++ii)
    {
        os <<  *ii << ",\n";
    }
    os << "]";
    return os;
}

template<typename T>
void getJsonPara( nlohmann::json j2, std::string para_name, T& defaultValue)
{
    try{
        if (j2.find(para_name) != j2.end()) {
            defaultValue = j2[para_name].get<T>();
            std::cout << "Read " << para_name << " = " << defaultValue << "\n";
        }
    }
    catch (std::exception& e){
        std::cerr << "exception reading parameter "<< para_name << "\n";
        std::cerr << e.what() << "\n";
        return ;
    }
}
