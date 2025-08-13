/*
La curva formada por uncable colgante se llama catenaria. Supongamos que el punto más bajo
de una catenaria es el origen (0; 0), entonces la ecuación de la catenaria es y = C * cosh(a/C) - C.
Si queremos determinar la catenaria que pasa por los puntos (±a; b), entonces debemos resolver la
ecuación b = C * cosh(a/C) - C donde la incógnita es C.

Inciso 1:
Demostrar que la catenaria que pasa por los puntos (±10, 6) es:
y = 9.1889 * cosh(a/9.1889) - 9.1889

X0 = 7, tolerance = 1e-4
El punto fijo de g(x) o la raiz de f(x) es igual a 9.1889
El error es igual a 0.0000
Resultado : 9.1889 +- 0.0000
Se han iterado 8 veces
Es decir que C = 9.1889, quedando demostrado el inciso 1.
Nos queda la catenaria: y = 9.1889 * cosh(a/9.1889) - 9.1889

Inciso 2:
Halle la catenaria que pasa por los puntos (±12; 5)
Debemos reemplazar en la ecuacion general de la catenaria, quedandonos 5 = C * cosh(12/C) - C, y despejando de aqui nos queda que
0 = C * cosh(12/C) - C - 5 = f(x) por lo tanto x = g(x) nos queda -> C = C * cosh(12/C) - 5

X0 = 16, tolerance = 1e-4
El punto fijo de g(x) o la raiz de f(x) es igual a 15.1672
El error es igual a 0.0001
Resultado : 15.1672 +- 0.0001
Se han iterado 19 veces
Es decir que C = 15.1672, quedandonos la catenaria igual a:
y = 15.1672 * cosh(a/15.1672) - 15.1672
Se verificó graficando esta última ecuación
*/

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

    printf("El punto fijo de g(x) o la raiz de f(x) es igual a %.4f\n", newX);
    printf("El error es igual a %.4f\n", error); 
    printf("Resultado : %.4f +- %.4f\n", newX, error);
    printf("Se han iterado %d veces\n", iteration);

    return 0;
}

double g(double x) {
    //return ((x*cosh(10.0/x)) - 6); Hecho para el primer inciso del ejercicio
    return ((x * cosh(12/x)) - 5);
}

double f(double x) {
    //return ((x*cosh(10/x)) - x - 6); Hecho para el primer inciso del ejercicio
    return ((x * cosh(12/x)) - x - 5); 
}