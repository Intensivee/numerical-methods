
#include <iostream>
#include <math.h>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip>
#include <limits.h>
#include <limits>

using namespace std;


vector<double> txtX;						// x pobrane z tabeli
vector<double> txtY;						// y pobrane z tabeli
vector<double> log10X;						// log10(x) pobrane z tabeli

vector<double> calculatedY;					// wartości funkcji obliczone za pomocą podstawowego wzoru
vector<double> betterCalculatedY;			// wartośći funkcji dla wzoru z wykorzystaniem szeregu taylora

vector<double> relativeError;				// |(W - Wp)/W| )		fla myFun
vector<double> betterRelativeError;			// |(W - Wp)/W| )		dla myFun2

vector<double> relativeErrorLog10;				// log10( |(W - Wp)/W| )		dla myFun
vector<double> betterRelativeErrorLog10;		// log10( |(W - Wp)/W| )		dla myFun2




// column == 0 wczytuje 1 kolumne, column == 1 wczytuje 2 kolumne, column == 2 wczytuje 3 kolumne
void readMyFIle(string fName, vector<double>* v, int column);

// funkcja zapisujaca vector do pliku
void saveData(string fName, vector<double>* v);



// f(x)
double myFun(double x) {
	return (1.0 - exp(-x)) / x;
}

// f(x) z wykorzystaniem szeregu taylora dla (1 - e^(-x)) (skróci sie odejmowanie ) a nastepnie caly uzyskany szereg skracamy na kartce z x (każdy el.)
double myFun2(double x ){
	//UWAGA LEPIEJ ZASTOSOWAC JEDNA ZMIENNA DO FACTORIAL I X POWER i w każdym kroku *=x/i;
	
	if (x < 1) {			
		double suma = 0.0;
		double factorial = 1.0;
		double xpower = 1;
		int sign = -1;

		for (int i = 1; i < 1000; i++) {
			factorial *= i;
			sign *= -1;
			suma += sign * xpower / factorial;
			xpower *= x;
		}
		return (double)(suma);
	}
	else return myFun(x);

}

//najlepsza metoda
double myFun3(double x) {
	if (x < 1) {
		double suma = 1.0;
		double previous = 1;
		int sign = 1;
		
		for (int i = 1; i < 1000; i++) {
			sign *= -1;
			previous *= x / (i+1);
			suma += sign * previous;
		}
		return (double) suma;
	}
	else return myFun(x);
}



int main()
{

	readMyFIle("data", &log10X, 0);
	readMyFIle("data", &txtX, 1);
	readMyFIle("data", &txtY, 2);
	
	// wyznaczam f(x)
	for (auto i : txtX) {
		calculatedY.push_back(myFun(i));
	}

	// lepsza metoda (z wykrozystaniem szeregu )
	int z = 0;
	for (auto i : txtX) {
		betterCalculatedY.push_back(myFun3(i));
		z++;
	}

	

	// obliczanie log10 z błędu względnego f(x)
	double temp;
	for (int i = 0; i < txtY.size(); i++) {
		temp = fabs((txtY[i] - calculatedY[i]) / txtY[i]);
		relativeError.push_back(temp);
		relativeErrorLog10.push_back(log10(temp));
	}

	for (int i = 0; i < txtY.size(); i++) {
		temp = fabs((txtY[i] - betterCalculatedY[i]) / txtY[i]);
		betterRelativeError.push_back(temp);
		betterRelativeErrorLog10.push_back(log10(temp));
	}

	
	
	saveData("log10X", &log10X);
	saveData("relativeErrorLog10", &relativeErrorLog10);
	saveData("relativeError", &relativeError);
	saveData("txtX", &txtX);
	saveData("txtY", &txtY);
	saveData("calculatedY", &calculatedY);


	
	// wypisuje X | Y | Y obliczony 2 metodą | blad wzgledny 2 metody
	for (int k = 0; k < log10X.size(); k++) {
		cout << fixed << setprecision(20) << scientific << txtX[k] << " ";
		cout << fixed << setprecision(20) << scientific << exp(-txtX[k]) << " ";
		cout << fixed << setprecision(20) << scientific << txtY[k] << " ";
		cout << fixed << setprecision(20) << scientific << betterCalculatedY[k] << " ";
		cout << fixed << setprecision(20) << scientific << betterRelativeError[k] << endl;
	}

	ofstream of1("wykres1.txt");		// wykres dla myFun()
	ofstream of2("wykres2.txt");		// wykres dla myFun2()  - wykorzystujacej szereg taylora
	for (int k = 0; k < log10X.size()-1; k++) {
		of1 << fixed << setprecision(20) << scientific << log10X[k] << " ";
		of1 << fixed << setprecision(20) << scientific << relativeErrorLog10[k] << endl;

		of2 << fixed << setprecision(20) << scientific << log10X[k] << " ";
		of2 << fixed << setprecision(20) << scientific << betterRelativeErrorLog10[k] << endl;
	}
	
}






void readMyFIle(string fName, vector<double>* v, int column) {
	try {
		string notImportant;
		double important;
		ifstream file("data.txt");
		if (file.fail()) throw runtime_error("blad otwarcia pliku");

		while (file.good()) {
			if (column == 0) {
				file >> important;
				file >> notImportant;
				file >> notImportant;
			}
			else if (column == 1) {
				file >> notImportant;
				file >> important;
				file >> notImportant;
			}
			else if (column == 2) {
				file >> notImportant;
				file >> notImportant;
				file >> important;
			}

			v->push_back(important);
		}
	}
	catch (runtime_error& e) {
		cout << e.what() << endl;
	}
	catch (...) {
		cout << "nieznany blad" << endl;
	}
}




void saveData(string fName, vector<double>* v) {
	ofstream of(fName + ".txt");
	for (auto i : *v) {
		of << fixed << setprecision(20) << scientific << i << endl;
	}
}


