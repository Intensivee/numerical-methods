#include <iostream>
#include <math.h>
#include <string>
#include <fstream>
#include <ostream>
#include <iomanip>

using namespace std;

string myFile = "data.txt";

// podstawowe funkcje
double fun_1(double x) {
	return sin(0.25 * x) * sin(0.25 * x) - x;
}
double fun_2(double x) {
	return tan(2 * x) - 1.0 - x;
}

// fi dla picarda
double fi_1(double x) {
	return sin(0.25 * x) * sin(0.25 * x);
}
double fi_2(double x) {
	return tan(2 * x) - 1.0;
}
// pochodna z fi dla picarda
double fi_prim_1(double x) {
	return 0.25 * sin(0.5 * x);
}
double fi_prim_2(double x) {
	return 2.0 / (cos(2.0*x) * cos(2.0*x));
}

// funkcje dla newtona ( f(x)/f'(x) )
double newton_fun1(double x) {
	return ((sin(0.25*x)*sin(0.25*x)) - x) / (0.25 * sin(0.5*x) - 1);
}
double newton_fun2(double x) {
	return (tan(2 * x) - x - 1) / (2 / (cos(2 * x)* cos(2 * x)) - 1);
}

// metoda picarda
void picard(string fname, double fi(double), double fi_prim(double), double f(double), double x_n, int n_max, double tol_x, double tol_f) {
	fstream o;
	// otwieram plik na dopisywanie lub tworzenie
	o.open(myFile, ios::out | ios::app);

	o << "------------------------------------  " << fname << "  --------------------------------------" << endl;
	o << fixed << setw(7) << "  n |" << setw(25) << "         x_n1           |" << setw(25) << "          e_n           |";
	o << setw(25) << "       residuum         |" << endl;

	double x_n1;
	double e_n;
	double residuum;

	// głowna pętla obliczająca kolejne przybliżenia (arbitralne ograniczenie na liczbę iteracji) 
	for (int n = 1; n <= n_max; n++) {

		// warunek zbieżności | fi'(x) | < 1 => zbieżność
		if (fabs(fi_prim(x_n)) < 1) {
			x_n1 = fi(x_n);
			e_n = fabs(x_n1 - x_n);
			x_n = x_n1;
			residuum = fabs(f(x_n));

			o << setprecision(16) << setw(5) << n << " |" << setw(23) << x_n << " |" << setw(23) << e_n << " |";
			o << setw(23) << residuum << " |" << endl;

			// kryterium dokładności wyznaczenia X*
			if (e_n <= tol_x) {
				o << "--------------------kryterium dokładności wyznaczenia Xn----------------- " << endl;
				break;
			}

			// kryterium wiarygoności Xn jako przybliżenia pierwiastka (residuum)
			if (residuum <= tol_f) {
				o << "--------------------kryterium wiarygodności Xn jako przybliżenia X*---------------- " << endl;
				break;
			}
		}
		else {
			o << "fi(x) nie jest zwezajace i nie zachodzi zbieznosc iteracji." << "\n\n\n\n" << endl;
			return;
		}
	}
	// ostatnie wyznaczone przybliżenie
	o << "f(x) = 0 ---> \nx = " << x_n << "\n\n\n\n";
	o.close();
}




void bisekcji(string fname, double f(double), double a, double b, int n_max, double tol_x, double tol_f) {
	fstream o;
	// otwieram plik na dopisywanie lub tworzenie
	o.open(myFile, ios::out | ios::app);

	o << "---------------------  " << fname << "  ---------------------------" << endl;
	o << fixed << setw(7) << "  n |" << setw(25) << "         x_n1           |" << setw(25) << "          e_n           |";
	o << setw(25) << "       residuum         |" << endl;

	//  warunek poprawności przedziału
	if (b < a) {
		o << " b < a - ERROR" << endl;
		return;
	}

	// warunek istnienia X* w CIĄGŁYM przedziale (a,b) 
	else if (f(b) * f(a) > 0) {
		o << "Funkcja na koncach przedzialu przyjmuje te same znaki" << endl;
		return;
	}

	// dodatkowy warunek sprawdzajcy czy przypadkiem nie trafiono pierwiastka 
	else if (f(b) == 0) {
		o << "miejsce zerowe b = " << b << endl;
		return;
	}
	else if (f(a) == 0) {
		o << "miejsce zerowe b = " << a << endl;
		return;
	}

	double x_n;
	double e_n;
	double residuum;

	// głowna pętla obliczająca kolejne przybliżenia (arbitralne ograniczenie na liczbę iteracji) 
	for (int n = 1; n <= n_max; n++) {

		x_n = (a + b) / 2;
		e_n = fabs((b - a) / 2);
		residuum = fabs(f(x_n));

		// sprawdzenie który przedział spełnia warunek miejsca zerowego 
		if (f(a) * f(x_n) < 0) b = x_n;
		else a = x_n;

		o << setprecision(16) << setw(5) << n << " |" << setw(23) << x_n << " |" << setw(23) << e_n << " |";
		o << setw(23) << residuum << " |" << endl;

		// kryterium dokładności wyznaczenia X*
		if (e_n <= tol_x) {
			o << "--------------------kryterium dokładności wyznaczenia Xn----------------- " << endl;
			break;
		}

		// kryterium wiarygoności Xn jako przybliżenia pierwiastka (residuum)
		if (residuum <= tol_f) {
			o << "--------------------kryterium wiarygodności Xn jako przybliżenia X*---------------- " << endl;
			break;
		}
	}
	// ostatnie wyznaczone przybliżenie
	o << "f(x) = 0 ---> x = " << x_n << "\n\n\n\n";
	o.close();
}




void newtona(string fname, double newton_fun(double), double f(double), double x_n, int n_max, double tol_x, double tol_f) {
	fstream o;
	// otwieram plik na dopisywanie lub tworzenie
	o.open(myFile, ios::out | ios::app);

	o << "---------------------  " << fname << "  ---------------------------" << endl;
	o << fixed << setw(7) << "  n |" << setw(25) << "         x_n            |" << setw(25) << "          e_n           |";
	o << setw(25) << "       residuum         |" << endl;

	double e_n;
	double residuum;
	double x_n1;

	// głowna pętla obliczająca kolejne przybliżenia (arbitralne ograniczenie na liczbę iteracji) 
	for (int n = 1; n <= n_max; n++) {

		x_n1 = x_n - newton_fun(x_n);
		residuum = fabs(f(x_n1));
		e_n = fabs(x_n1 - x_n);
		x_n = x_n1;

		o << setprecision(16) << setw(5) << n << " |" << setw(23) << x_n1 << " |" << setw(23) << e_n << " |";
		o << setw(23) << residuum << " |" << endl;

		// kryterium dokładności wyznaczenia X*
		if (e_n <= tol_x) {
			o << "--------------------kryterium dokładności wyznaczenia Xn----------------- " << endl;
			break;
		}

		// kryterium wiarygoności Xn jako przybliżenia pierwiastka (residuum)
		if (residuum <= tol_f) {
			o << "--------------------kryterium wiarygodności Xn jako przybliżenia X*---------------- " << endl;
			break;
		}
	}
	// ostatnie wyznaczone przybliżenie
	o << "f(x) = 0 ---> x = " << x_n << "\n\n\n\n";
	o.close();
}





void siecznych(string fname, double f(double), double x_n, double x_n1, int n_max, double tol_x, double tol_f) {
	fstream o;
	// otwieram plik na dopisywanie lub tworzenie
	o.open(myFile, ios::out | ios::app);

	o << "---------------------  " << fname << "  ---------------------------" << endl;
	o << fixed << setw(7) << "  n |" << setw(25) << "         x_n2           |" << setw(25) << "          e_n           |";
	o << setw(25) << "       residuum         |" << endl;

	double e_n;
	double residuum;
	double x_n2;
	// głowna pętla obliczająca kolejne przybliżenia (arbitralne ograniczenie na liczbę iteracji) 
	for (int n = 1; n <= n_max; n++) {

		x_n2 = x_n1 - (f(x_n1) / ((f(x_n1) - f(x_n)) / (x_n1 - x_n)));
		residuum = fabs(f(x_n2));
		e_n = fabs(x_n2 - x_n1);
		x_n = x_n1;
		x_n1 = x_n2;

		o << setprecision(16) << setw(5) << n << " |" << setw(23) << x_n2 << " |" << setw(23) << e_n << " |";
		o << setw(23) << residuum << " |" << endl;

		// kryterium dokładności wyznaczenia X*
		if (e_n <= tol_x) {
			o << "--------------------kryterium dokładności wyznaczenia Xn----------------- " << endl;
			break;
		}

		// kryterium wiarygoności Xn jako przybliżenia pierwiastka (residuum)
		if (residuum <= tol_f) {
			o << "--------------------kryterium wiarygodności Xn jako przybliżenia X*---------------- " << endl;
			break;
		}
	}
	// ostatnie wyznaczone przybliżenie
	o << "f(x) = 0 ---> x = " << x_n2 << "\n\n\n\n";
	o.close();
}

int main()
{

	picard("Metoda_Picarda_Fun1, xn=6", fi_1, fi_prim_1, fun_1, 6, 100, 0.000000000000001, 0.000000000000001);
	picard("Metoda_Picarda_Fun2, xn=0.6", fi_2, fi_prim_2, fun_2, 0.6, 10, 0.00000000001, 0.00001);
	bisekcji("Metoda_Bisekcji_Fun1, a=-10, b=111", fun_1, -10, 111, 100, 0.000000000000001, 0.000000000000001);
	bisekcji("Metoda_Bisekcji_Fun2, a=-0.4, b=0.5", fun_2, -0.4, 0.5, 1000, 0.000000000000001, 0.000000000000001);
	newtona("Metoda_Newtona_Fun1, xn=1000", newton_fun1, fun_1, 1000, 100, 0.000000000000001, 0.000000000000001);
	newtona("Metoda_Newtona_Fun2, xn=0.6", newton_fun2, fun_2, 0.6, 100, 0.000000000000001, 0.000000000000001);
	siecznych("Metoda_Siecznych_Fun1, xn=0.6, xn1=0.4", fun_1, 0.6, 0.4, 100, 0.000000000000001, 0.000000000000001);
	siecznych("Metoda_Siecznych_Fun2, xn=0.6, xn1=0.5", fun_2, 0.6, 0.5, 100, 0.000000000000001, 0.000000000000001);
}

