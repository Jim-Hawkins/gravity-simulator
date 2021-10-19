#include <iostream>
#include <limits>
#include <iomanip>

using namespace std;

int main(){
	double m = 9.9999993774659967e+20;
	double F = 3.7111137423073675e+24;
	cout << std::setprecision(std::numeric_limits<double>::max_digits10) << F << endl;
	cout << std::setprecision(std::numeric_limits<double>::max_digits10) << m << endl;
	cout << std::setprecision(std::numeric_limits<double>::max_digits10)<< F/m << endl;
	cout << std::setprecision(std::numeric_limits<double>::max_digits10)<< 1/m*F << endl;
	return 0;
}
