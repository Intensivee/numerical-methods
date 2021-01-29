#define _USE_MATH_DEFINES
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
using namespace std;

double fun(double x) {
	return 1.0 / (1.0 + 10.0*x*x*x*x*x*x);
}
double newton(double x, double* Xi, double* Fi);
void getRownolegle(double* Xi, double* Fi);
void getCzybyszew(double* Xi, double* Fi);
double* get_c(double * X, double *F);


// przedział [a;b]
static double a = -1.0;
static double b = 1.0;

// ilosc węzłów
static int k = 22;
// ilość iteracji po przedziale [-1;1]
static int N = 100;

int main()
{
	
	// ogarnij inny krok dla funkcji niż węzły? xd
	ofstream o("data.txt");
	double* X_czybyszew = new double[k];
	double* F_czybyszew = new double[k];
	double* X_rownolegle = new double[k];
	double* F_rownolegle = new double[k];
	double* X_1 = new double[k];
	double* X_2 = new double[k];
	double* c_rownolegle = new double[k];
	double* c_czybyszew = new double[k];



	double h = (b - a) / (N - 1);

	getRownolegle(X_rownolegle, F_rownolegle);
	getCzybyszew(X_czybyszew, F_czybyszew);
	c_rownolegle = get_c(X_rownolegle, F_rownolegle);
	c_czybyszew = get_c(X_czybyszew, F_czybyszew);


	double x = a;
	for (int i = 0; i < N; i++, x += h) {
		copy(X_rownolegle, X_rownolegle + k, X_1);
		copy(X_czybyszew, X_czybyszew + k, X_2);

		o << fixed << setprecision(16) << x << " " << fun(x) << " " << newton(x, X_1, c_rownolegle) << " " << newton(x, X_2, c_czybyszew) << endl;
	}

}



double* get_c(double * X, double *F) {
	double* c = new double[k];
	// wyznaczam wsp. c. r - rząd ilorazu
	c[0] = F[0];
	for (int r = 1; r < k; r++) {

		// wyznaczam kolejne kolumny tabeli ilorazów róznicowych rzędu r
		for (int j = 0; j < k - r; j++) {
			F[j] = (F[j + 1] - F[j]) / (X[r + j] - X[j]);
		}

		// zapisuje wsp. c (wartość na górze kolumny)
		c[r] = F[0];
	}
	return c;
}

double newton(double x, double* X, double* c) {

	double val = c[k-1];
	double* polynomials = new double[k];

	// wyznaczam kolejne wsp. wielomianu dla arg. x
	for (int i = 0; i < k; i++) {
		polynomials[i] = 1.0;

		for (int j = 0; j < i; j++) {
			polynomials[i] *= x - X[j];
		}
	}

	// obliczam wartość funkcji (horner)
	for (int i = k-1; i > 0; i--) {
		val = val * (x - X[i - 1]) + c[i - 1];
	}

	delete polynomials;
	return val;
}


void getCzybyszew(double* X, double* F) {
	for (int i = 0; i < k; i++) {
		X[i] = (b + a) / 2.0 + (b - a) / 2.0 * cos((2.0 * i + 1.0) / (2.0 * (k - 1.0) + 2.0) * M_PI);
		F[i] = fun(X[i]);
	}
}

void getRownolegle(double* X, double* F) {
	double val = a;

	// krok co jaki ustawiam węzeł
	double h2 = (b - a) / (k - 1.0);
	for (int i = 0; i < k; i++, val += h2) {
		cout << "i: " << i << "\t x: " << val << endl;
		X[i] = val;
		F[i] = fun(val);
	}
}