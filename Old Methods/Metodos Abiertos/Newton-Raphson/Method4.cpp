#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

using namespace std;

double fPrima(double x); //Es la derivada de f(x)
double f(double x); 

int main(int argc, char const *argv[]) {
    double tolerance;
    double error = 1;
    double X0; //Generalmente es 0
    double oldX;
    double newX;
    int iteration = 0;

    cout<<"Ingrese el punto del eje x desde donde empezará a buscar la raíz"<<endl;
    cin>>X0;
    cout<<"Ingrese la tolerancia deseada"<<endl;
    cin>>tolerance;

    oldX = X0;
    do {
        iteration++;
        if(fabs(fPrima(oldX)) < 1e-5) {
            cout<<"Derivada muy pequeña"<<endl;
            exit(0);
        } else {
            newX = oldX - (f(oldX) / fPrima(oldX));
            error = fabs(newX - oldX);
            oldX = newX;
        }
    } while((error > tolerance) && (iteration < 10000));

    printf("La raiz de f(x) es igual a %.5f\n", newX);
    printf("El error es igual a %.5f\n", error); 
    printf("Resultado : %.5f +- %.5f\n", newX, error);
    printf("Se han iterado %d veces\n", iteration);

    return 0;
}

double fPrima(double x) {
    return (pow(x, 9) * 10);
}

double f(double x) {
    return (pow(x, 10) - 1);
}