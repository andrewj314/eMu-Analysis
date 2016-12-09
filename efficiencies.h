#ifndef efficiencies_h
#define efficiencies_h
#include <map>
#include <utility>
#include <string>

template <typename T1, typename T2> std::pair<double, std::pair<double, double> > GetEfficiency(T1, T2);
template <typename T1, typename T2> std::string PrintEfficiency(T1, T2);
std::pair<double, double> GetSurvivingEvents(std::pair<double, double> nEvents, 
                          std::pair<double, std::pair<double, double> > efficiency);
#endif // #ifdef 
