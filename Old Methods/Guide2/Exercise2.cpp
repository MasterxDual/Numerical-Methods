/*
Realice un programa adecuado para evaluar el punto fijo de funciones arbitrarias. Utilice este
programa para aproximar los puntos fijos (si es que hay alguno) de cada una de las siguientes
funciones. Las respuestas deben tener 12 cifras decimales exactas. Grafique además cada función y
la recta y = x para poder visualizar claramente los puntos fijos si es que existen.

g(x) = pow(x, 5) - 3*pow(x, 3) - 2*pow(x, 2) + 2
Nunca converge. Su derivada es mayor a 1.

g(a) = cos(sin(a)), f(a) = cos(sin(a));
El punto fijo de g(a) o la raiz de f(a) es igual a 0.768189133355
El error es igual a 0.000063361005
Resultado : 0.768189133355 +- 0.000063361005
Se han iterado 12 veces

g(n) = pow(x, (x - cos(x))), f(n) = pow(x, (x - cos(x)));
X0 = 1.1, tolerance = 1e-4
El punto fijo de g(n) o la raiz de f(n) es igual a 1.000043128758
El error es igual a 0.000050658209
Resultado : 1.000043128758 +- 0.000050658209
Se han iterado 11 veces
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

    printf("El punto fijo de g(n) o la raiz de f(n) es igual a %.12f\n", newX);
    printf("El error es igual a %.12f\n", error); 
    printf("Resultado : %.12f +- %.12f\n", newX, error);
    printf("Se han iterado %d veces\n", iteration);

    return 0;
}

double g(double x) {
   //return (pow(x, 5) - (3*pow(x, 3)) - (2*pow(x, 2)) + 2); Hemos utilizado la g(x) que ya viene por defecto de la guia

   //return (cos(sin(x))); Hemos utilizado la g(a) que ya viene por defecto de la guia
   
   return (pow(x, (x - cos(x)))); //Hemos utilizado la g(n) que ya viene por defecto de la guia
}

double f(double x) {
    //return (pow(x, 5) - 3*pow(x, 3) - 2*pow(x, 2) + 2); g(x) que viene por defecto de la guia, la utilice como f(x) y a su vez como g(x) ya que no especifica en la guia de ejercicios

    //return (cos(sin(x))); g(a) que viene por defecto de la guia, la utilice como f(a) y a su vez como g(a) ya que no especifica en la guia de ejercicios

    return (pow(x, (x - cos(x)))); //g(n) que viene por defecto en la guia, la utilice como f(n) y a su vez como g(n) ya que no se especifica en la guia de ejercicios
}