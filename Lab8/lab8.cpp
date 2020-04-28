#define _USE_MATH_DEFINES
#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>

#define N 30

using namespace std;

// dzięki skorzystaniu z osobnej funkcji program jest otwarty na modyfikacje
template <class T> T f(T x) {
	return sin(x);
}

// wartość dokładna 
template <class T> T f_prim(T x) {
	return cos(x);
}

// różniczki
template <class T> T forward_difference(T x, T h) {
	return (f(x + h) - f(x)) / h;
}

template <class T> T backward_difference(T x, T h) {
	return (f(x) - f(x - h)) / h;
}

template <class T> T central_difference(T x, T h) {
	return (f(x + h) - f(x - h)) / (T(2) * h);
}

// wyciągam 1/2 przed nawias w celu zmniejszenia ilości operacji 
template <class T> T three_point_forward_difference(T x, T h) {
	return (T(-3) * f(x) + T(4) * f(x + h) - f(x + T(2) * h)) / (T(2) * h);
}

// wyciągam 1/2 przed nawias w celu zmniejszenia ilości operacji 
template <class T> T three_point_backward_difference(T x, T h) {
	return (f(x - T(2) * h) - T(4) * f(x - h) + T(3) * f(x)) / ( T(2) * h);
}

// liczby rzutuje na typ T ( tak aby odpowiedno były liczbą double lub float)
template <class T> void count_differences(string fname) {
	ofstream o(fname);
	o << fixed << setprecision(16);
	T x1 = T(0);
	T x2 = M_PI / T(4);
	T x3 = M_PI / T(2);
	T step = T(1);


	// log10( blad bezwzględny )
	for (int i = 1; i <= N; i++) {
		step /= T(10);
		o << log10(step) << " ";

		// różniczki na początku przedziału (x1)
		o << log10(fabs(f_prim(x1) - three_point_forward_difference(x1, step))) << " ";
		o << log10(fabs(f_prim(x1) - forward_difference(x1, step))) << " ";;

		// różniczki w środku przedziału (x2)
		o << log10(fabs(f_prim(x2) - backward_difference(x2, step))) << " ";
		o << log10(fabs(f_prim(x2) - central_difference(x2, step))) << " ";
		o << log10(fabs(f_prim(x2) - forward_difference(x2, step))) << " ";

		// różniczki na końcu przedziału (x3)
		o << log10(fabs(f_prim(x3) - backward_difference(x3, step))) << " ";
		o << log10(fabs(f_prim(x3) - three_point_backward_difference(x3, step))) << endl;
		
	}
	o.close();
}



int main()
{
	count_differences<double>("double_data.txt");
	count_differences<float>("float_data.txt");
}

