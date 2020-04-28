#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#define N 100000
using namespace std;

// ponieważ macierz trój-diagonalna jest macierzą rzadką, zamiast tworzyć tablie kwadratową korzystam z wektorów
void Thomas_D_reduction(double*L, double* D, double* U, int n);
void Thomas_solve(double* L, double* D, double* U, double* b, double* x, int n);
void get_convential_task_vectors(double* L, double* D, double* U, double* b, double x, double h, int n);
void get_numerow_task_vectors(double* L, double* D, double* U, double* b, double x, double h, int n);
void discretization(string fname, void method(double*, double*, double*, double*, double, double, int));
double vector_max(double* v, int n);

// przedział x
const double x_start = 0.0;
const double x_end = 1.0;

// wsp. p,q,r,s
const double p = 1.0;
const double q = 0.0;
const double r = -4.0;
double s_fun(double x) {
	return -x;
}

// wssp. warunku brzegowego U(0)=1
const double alfa = 0.0;
const double beta = 1.0;
const double gamma = -1.0;

// wsp. warunku brzegowego U(1)=0
const double phi = 0.0;
const double psi = 1.0;
const double theta = 0.0;

double analitical_value(double x) {
	return (exp(2.0 - 2.0*x) - 4.0*exp(4.0 - 2.0*x) + 4.0*exp(2.0*x) - exp(2.0 + 2.0*x) - x + x * exp(4.0)) / (4.0 - 4.0*exp(4.0));
}

 
int main()
{
	discretization("konwencjonalna", get_convential_task_vectors);
	discretization("Numerowa", get_numerow_task_vectors);
}


void discretization(string fname, void method(double*, double*, double*, double*, double, double, int)) {

	ofstream o1("dyskretyzacja_" + fname + ".txt");
	ofstream o2("log10_blad_dyskretyzacja_" + fname + ".txt");

	// przybliżone wartości, blad
	double* counted_x = new double[N];
	double* error = new double[N];

	// alokuje pamięć jedynie raz 
	double* L = new double[N];
	double* D = new double[N];
	double* U = new double[N];
	double* b = new double[N];
	double h;
	double x_temp;

	// ilość kroków  
	for (int n = 50; n <= N; n+=50) {

		// obliczam krok h
		h = (x_end - x_start) / (n-1);

		// Macierz A przedstawiam za pomocą wektorów L,D,U
		method(L, D, U, b, x_start, h, n);

		// alg. Thomasa
		Thomas_D_reduction(L, D, U, n);
		Thomas_solve(L, D, U, b, counted_x, n);

		// wyznaczam blad przyblizenia
		x_temp = x_start;
		for (int i = 0; i < n; i++) {
			error[i] = fabs(counted_x[i] - analitical_value(x_temp));
			x_temp += h;
		}

		// znajduje max blad i wpisuje go do pliku
		o2 << setprecision(16) << log10(h) << " " << log10(vector_max(error, n)) << endl;

		// jednorazowo wyznaczam wartości analityczne i przybliżone do wykresu
		// n==100 <-- większa ilość punktów powodowałaby nieczytelność takiego wykresu
		if (n == 100) {
			x_temp = x_start;
			for (int i = 0; i < n; i++) {
				o1 << setprecision(16) << x_temp << " " << counted_x[i] << " " << analitical_value(x_temp) << endl;
				x_temp += h;
			}
		}
	}
	o1.close();
	o2.close();
	delete[] counted_x;
	delete[] error;
	delete[] L;
	delete[] D;
	delete[] U;
	delete[] b;
}

void get_convential_task_vectors(double* L, double* D, double* U, double* b, double x, double h, int n) {

	L[0] = 0.0;
	D[0] = beta - alfa / h;
	U[0] = alfa / h;
	b[0] = -gamma;

	for (int i = 1; i < n - 1; i++) {
		L[i] = p / (h*h) - q / (2.0*h);
		D[i] = r - (2.0 * p) / (h*h);
		U[i] = p / (h*h) + q / (2.0*h);
		b[i] = -s_fun(x + i * h);
	}

	L[n - 1] = -phi / h;
	D[n - 1] = phi / h + psi;
	U[n - 1] = 0.0;
	b[n - 1] = -theta;
}

void get_numerow_task_vectors(double* L, double* D, double* U, double* b, double x, double h, int n) {

	L[0] = 0.0;
	D[0] = beta - alfa / h;
	U[0] = alfa / h;
	b[0] = -gamma;

	for (int i = 1; i < n - 1; i++) {
		L[i] = p / (h * h) + r / 12.0;
		D[i] = (-2.0 * p) / (h * h) + r * (10.0 / 12.0);
		U[i] = p / (h * h) + r / 12.0;
		b[i] = -s_fun(x + i * h - h) / 12.0 - s_fun(x + i * h)*10.0 / 12.0 - s_fun(x + i * h + h) / 12.0;
	}

	L[n - 1] = -phi / h;
	D[n - 1] = phi / h + psi;
	U[n - 1] = 0.0;
	b[n - 1] = -theta;
}

void Thomas_D_reduction(double*L, double* D, double* U, int n) {
	// D[0] pozostawiam bez zmian
	for (int i = 1; i < n; i++)
		D[i] -= L[i] / D[i - 1] * U[i - 1];
}

void Thomas_solve(double* L, double* D, double* U, double* b, double* x, int n) {
	// redukuje macierz wyników b->r
	// b[0] pozostawiam bez zmian
	for (int i = 1; i < n; i++)
		b[i] -= L[i] / D[i - 1] * b[i - 1];

	// wyznaczam x
	x[n - 1] = b[n - 1] / D[n - 1];
	for (int i = N - 2; i >= 0; i--)
		x[i] = (b[i] - U[i] * x[i + 1]) / D[i];
}


double vector_max(double* v, int n) {
	double max = v[0];
	for (int i = 1; i < n; i++) {
		if (v[i] > max)
			max = v[i];
	}
	return max;
}

