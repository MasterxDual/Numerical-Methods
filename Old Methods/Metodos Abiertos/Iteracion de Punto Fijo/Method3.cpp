#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

using namespace std;

double g(double x); //A g(x) la conseguimos despejando x de f(x)
double f(double x); //No la utilizamos, solo la dejo plasmada para mostrar la funcion desde la que partimos

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
        if(fabs((g(oldX+0.01) - g(oldX)) / 0.01) >= 1) {
            cout<<"No se cumple el criterio de convergencia"<<endl;
            exit(0);
        } else {
            newX = g(oldX);
            error = fabs(oldX - newX);
            oldX = newX;
        }
    } while(error > tolerance);

    printf("El punto fijo de g(x) o la raiz de f(x) es igual a %.5f\n", newX);
    printf("El error es igual a %.5f\n", error); 
    printf("Resultado : %.5f +- %.5f\n", newX, error);
    printf("Se han iterado %d veces\n", iteration);

    return 0;
}

double g(double x) {
    return ((-1/4) * sin(x));
}

double f(double x) {
    return ((1/2) * sin(x) + (2*x));
}