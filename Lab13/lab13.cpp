#define _USE_MATH_DEFINES
#include <iostream>
#include <math.h>
#include <iomanip>
#include "calerf.h"

using namespace std;

double fun(double x) {
	return calerf::ERFCL(x);
}

double inside_fun(double y) {
	return exp(-(y*y));
}

double rectangular_left(double a, double b, double dx);
double rectangular_right(double a, double b, double dx);
double rectangular_midle(double a, double b, double dx);
double trapezoidal(double a, double b, double dx);
double parabol(double a, double b, double dx);

int main()
{
	cout << fixed << setprecision(16);
	double a = -1.0;
	double b = 1.0;
	int N = 100;
	double dx = (b - a) / N;

	double pre = 2.0 / pow(M_PI, 0.5);
	cout << fun(2);
	cout << pre * rectangular_left(a, b, dx) << endl;
	cout << pre*rectangular_right(a, b, dx) << endl;
	cout << pre*rectangular_midle(a, b, dx) << endl;
	cout << pre*trapezoidal(a, b, dx) << endl;
	cout << pre*parabol(a, b, dx) << endl;
}



// metoda prostokątów wariant z węzłem interpolacji po lewej stronie przedziału
double rectangular_left(double a, double b, double dx)
{
	double integral = 0.0;
	// głowna pętla wyznaczająca sume kolejnych f(a) (lewe wyrazy)
	for (double x = a; x < b; x += dx)
		integral += inside_fun(x);

	// wyciągam dx przed nawias
	integral *= dx;
	return integral;
}


// metoda prostokątów wariant z węzłem interpolacji po prawej stronie przedziału
double rectangular_right(double a, double b, double dx) {

	double integral = 0.0;
	// głowna pętla wyznaczająca sume kolejnych f(b) (prawe wyrazy)
	for (double x = a + dx; x <= b; x += dx)
		integral += inside_fun(x);

	// wyciągam dx przed nawias
	integral *= dx;
	return integral;
}


// metoda prostokątów wariant z węzłem interpolacji po środku przedziału
double rectangular_midle(double a, double b, double dx)
{
	double integral = 0.0;
	// głowna pętla wyznaczająca sume kolejnych f((a+b)/2)
	for (double x = a + dx/2.0; x < b; x += dx)
		integral += inside_fun(x);

	// wyciągam dx przed nawias
	integral *= dx;
	return integral;
}


// metoda trapezów
double trapezoidal(double a, double b, double dx) {

	double integral = 0.0;
	// głowna pętla wyznaczająca sume kolejnych wyrazów poza pierwszym i ostatnim
	// (wszystkie wyrazy poza 1 i ostatnim występują 2 krotnie)
	for (double x = a + dx; x < b; x += dx)
		integral += inside_fun(x);

	// pojawiają się tylko raz więc dzielenie się nie skraca
	integral += (inside_fun(a) + inside_fun(b)) / 2.0;
	integral *= dx;
	return integral;
}


// metoda parabol
double parabol(double a, double b, double dx) {

	double middle = 0.0;
	double edge = 0.0;
	// głowna pętla wyznaczająca sume kolejnych wyrazów poza pierwszym i ostatnim
	// (wszystkie wyrazy poza pierwszym, ostatnim i wszystkimi (x+dx/2 - środkowymi) występują 2 krotnie)
	for (double x = a + dx; x <= b; x += dx) {
		middle += inside_fun(x - dx / 2);

		if(x < b)
			edge += inside_fun(x);
	}

	double integral = (inside_fun(a) + inside_fun(b) + 2.0 * edge + 4.0 * middle) / 6.0 * dx;
	return integral;
}