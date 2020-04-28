#include <iostream>
#include <math.h>
#include <string>
#include <fstream>
#include <ostream>
#include <iomanip>

using namespace std;


double my_fun_1(double x, double y, double z) {
	return x * x + y * y + z * z - 2.0;
}

double my_fun_2(double x, double y) {
	return x * x + y * y - 1;
}

double my_fun_3(double x, double y) {
	return x * x - y;
}

// Dzięki funkcyjnemu zwracaniu Odwróconej macierzy Jacobianu program łatwiej przerobić na inne równania
double** Inverse_Jacobi_Matrix(double** J, double x, double y, double z) {

	// odwrotnosc wyznacznika macierzy Jakobiego
	double inverse_det = 1.0 / (-4.0 * x * z - 8.0 * x * y * z);

	// Wyznaczenie macierzy odwrotnej do Jacobiego, dla danych x,y,z 
	J[0][0] = 0;
	J[0][1] = inverse_det * -2.0 * z ;
	J[0][2] = inverse_det * -4.0 * y * z ;
	J[1][0] = 0;
	J[1][1] = inverse_det * -4.0 * x * z;
	J[1][2] = inverse_det * 4.0 * x * z;
	J[2][0] = inverse_det * (-2.0 * x - 4.0 * x * y);
	J[2][1] = inverse_det * (2.0 * x + 4.0 * x * y);
	J[2][2] = 0;
	return J;
	
}

double* calculate_fun_values(double *V, double x, double y, double z) {
	V[0] = my_fun_1(x, y, z);
	V[1] = my_fun_2(x, y);
	V[2] = my_fun_3(x, y);
	return V;
}


double** alloc_matrix(int n, int m) {
	double **M = NULL;
	M = new double*[n];
	for (int i = 0; i < 3; i++)
		M[i] = new double[m];
	return M;
}

double* alloc_vector(int n) {
	double *V = NULL;
	V = new double[n];
	for (int i = 0; i < n; i++)
		V[i] = 0;
	return V;
}

double max_of_3(double a, double b, double c) {
	double max = (a < b) ? b : a;
	return ((max < c) ? c : max);
}

void Newton(string fname, double x, double y, double z, int n_max, double tol_x, double tol_f) {
	double e_n;
	double residuum;
	ofstream o(fname);

	// stworzenie odwróconej macierzy jacobiego, wektora wartości funkcji, oraz wektora kolejnego przybliżenia x,y,z
	double** inv_Jac = alloc_matrix(3, 3);
	double* fun_values = alloc_vector(3);
	double* X_n1 = alloc_vector(3);

	// obliczam wartości funkcji
	fun_values = calculate_fun_values(fun_values, x, y, z);



	o << "------------------------------------  " << fname << "  --------------------------------------" << endl;
	o << fixed << setw(7) << "  n |" << setw(25) << "          x             |" << setw(25) << "          y             |";
	o << setw(25) << "          z             |" << setw(25) << "          e_n           |";
	o << setw(25) << "       residuum         |" << endl;
	// Główna pętla licząca kolejne przybliżenia
	for (int i = 1; i <= n_max; i++) {

		// sprawdzam czy nie wystąpi dzielenie przez 0
		if (x == 0 || z == 0) {
			cout << "error" << endl;
			o << "--------------------dzielenie przez 0----------------- " << endl;
			// omijam wypisanie informacji o zakończeniu z powodu przekroczenia ilości iteracji
			goto END;
		}

		// uzupełnienie macierzy jacobiego 
		inv_Jac = Inverse_Jacobi_Matrix(inv_Jac, x, y, z);


		// Wyznaczanie kolejnego przybliżenia pierwiastka
		X_n1[0] = x;
		X_n1[1] = y;
		X_n1[2] = z;
		for (int n = 0; n < 3; n++) {
			for (int m = 0; m < 3; m++) {
				X_n1[n] = X_n1[n] - (inv_Jac[n][m] * fun_values[m]);
			}
		}

		// wyznaczam estymator błedu
		e_n = max_of_3(fabs(X_n1[0] - x), fabs(X_n1[1] - y), fabs(X_n1[2] - z));

		// przypisuje nowe wartości x,y,z
		x = X_n1[0];
		y = X_n1[1];
		z = X_n1[2];

		// obliczam wartości funkcji
		fun_values = calculate_fun_values(fun_values, x, y, z);

		// wyznaczam residuum
		residuum = max_of_3(fabs(fun_values[0]), fabs(fun_values[1]), fabs(fun_values[2]));


		// zapisywanie do pliku
		o << setprecision(16) << setw(5) << i << " |" << setw(23) << x << " |" << setw(23) << y << " |";
		o << setw(23) << z << " |" << setw(23) << e_n << " |";
		o << setw(23) << residuum << " |" << endl;


		// warunek tol_x
		if (e_n < tol_x) {
			o << "--------------------Spełniono Warunek dokladnosci X*----------------- " << endl;
			o << "f(X) = 0 ---> \nx = " << x << "\ny = " << y <<"\nz = "<< z << "\n\n\n\n";
			o.close();
			break;
		}
		// warunek residuum
		else if (residuum < tol_f) {
			o << "--------------------Spełniono Warunek residuum----------------- " << endl;
			o << "f(X) = 0 ---> \nx = " << x << "\ny = " << y << "\nz = " << z << "\n\n\n\n";
			o.close();
			break;
		}
	}
	
	// warunek ilosci iteracji
	o << "--------------------Warunek ilosci iteracji----------------- " << endl;
	o << "f(X) = 0 ---> \nx = " << x << "\ny = " << y << "\nz = " << z << "\n\n\n\n";

	END:
	o.close();
}

int main()
{
	Newton("data.txt",1.0, 1.0, 1.0, 100, 0.000000000000001, 0.000000000000001);
}

