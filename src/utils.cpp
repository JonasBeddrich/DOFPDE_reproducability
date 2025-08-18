#include "mfem.hpp"

using namespace std;
using namespace mfem;

void ReadWeightsPoles(std::vector<real_t> &weights, std::vector<real_t> &poles, const char *filename){
    std::ifstream file(filename);
    std::string line, item; 
    
    if (std::getline(file, line)) {
       std::stringstream ss(line);
       while (std::getline(ss, item, ',')) {
          weights.push_back(std::stod(item));
       }
    }
 
    if (std::getline(file, line)) {
       std::stringstream ss(line);
       while (std::getline(ss, item, ',')) {
          poles.push_back(-1.0 * std::stod(item));
       }
    }
 } 