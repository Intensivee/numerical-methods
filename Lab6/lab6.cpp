#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <sstream>
#define N 6

using namespace std;
// ponieważ macierz trój-diagonalna jest macierzą rzadką, zamiast tworzyć tablie kwadratową korzystam z wektorów
void set_task_vectors(double* L, double* D, double* U, double* b);
void Thomas_D_reduction(double*L, double* D, double* U);
void Thomas_solve(double* L, double* D, double* U, double* b, double* x);
void print_vector(double *V, string what);

int main()
{
	double* L = new double[N];
    double* D = new double[N];
    double* U = new double[N];
    double* b = new double[N];
    double* x = new double[N];
	set_task_vectors(L, D, U, b);
	Thomas_D_reduction(L, D, U);
	Thomas_solve(L, D, U, b, x);
	print_vector(x,"x");
}

// w celu zwiększenia mobilności programu (mail) dane wpisuje kolejno, zamiast zczytywać w pętli z pliku
void set_task_vectors(double* L, double* D, double* U, double* b) {
	L[0] = 0.0;
	L[1] = 1.0 / 3.0;
	L[2] = 1.0 / 5.0;
	L[3] = 1.0 / 7.0;
	L[4] = 1.0 / 9.0;
	L[5] = 1.0 / 11.0;

	D[0] = 10.0;
	D[1] = 20.0;
	D[2] = 30.0;
	D[3] = 30.0;
	D[4] = 20.0;
	D[5] = 10.0;

	U[0] = 1.0 / 2.0;
	U[1] = 1.0 / 4.0;
	U[2] = 1.0 / 6.0;
	U[3] = 1.0 / 8.0;
	U[4] = 1.0 / 10.0;
	U[5] = 0;

	b[0] = 31.0;
	b[1] = 165.0 / 4.0;
	b[2] = 917.0 / 30.0;
	b[3] = 851.0 / 28.0;
	b[4] = 3637.0 / 90.0;
	b[5] = 332.0 / 11.0;
}

void Thomas_D_reduction(double*L, double* D, double* U) {
	// D[0] pozostawiam bez zmian
	for (int i = 1; i < N; i++) 
		D[i] -= L[i] / D[i - 1] * U[i - 1];
}


void Thomas_solve(double* L, double* D, double* U, double* b, double* x) {
	// redukuje macierz wyników b->r
	// b[0] pozostawiam bez zmian
	for (int i = 1; i < N; i++)
		b[i] -= L[i] / D[i - 1] * b[i - 1];

	// wyznaczam x
	x[N - 1] = b[N - 1] / D[N - 1];
	for (int i = N - 2; i >= 0; i--)
		x[i] = (b[i] - U[i] * x[i + 1]) / D[i];
}


void print_vector(double *V, string what) {
	fstream o;
	o.open("data.txt", ios::out | ios::app);
	o << "------------------- " << what << " -------------------------" << endl;

	for (int i = 0; i < N; i++)
		o << V[i] << endl;
	o << "\n\n\n";
}