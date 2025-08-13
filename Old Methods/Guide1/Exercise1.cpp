/*Determine las raíces reales de:
f (x) = −2 + 7x − 5x2 + 6x3 
1. Un método de utilidad y que muchas veces sirve como guía para aproximar el valor de las
raíces de una ecuación determinada, es graficar la función para visualizar las raíces. Al graficar
la función se logra acotar el dominio de búsqueda. Utilizando algún graficador de su agrado,
obtenga una estimación de las raíces de f (x) de manera gráfica.

1. Graficamente podemos decir que la raiz se encuentra en el punto x = 0.330354

2. Utilizando el método de bisección para encontrar la raíz más pequeña de f(x). Adopte como 
valores iniciales xl = 0 y xu = 1, realice el proceso iterativo hasta que el error de la
aproximación se encuentre por debajo de 1 × 10−4.

2. a = 0, b = 1, tolerance = 1e-4;
La raiz es 0.33337
El error es 0.00006
Raiz: 0.33337 +- 0.00006
Se ha iterado 13 vez/veces
*/

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
using namespace std;

double f(double x);

int main(int argc, char const *argv[]) {
    double a;
    double b;
    double tolerance;
    double root;

    cout<<"Ingrese el limite izquierdo del intervalo de donde esta la raiz"<<endl;
    cin>>a;
    cout<<"Ingrese el limite derecho del intervalo de donde esta la raiz"<<endl;
    cin>>b;
    cout<<"Ingrese la tolerancia deseada del error"<<endl;
    cin>>tolerance;

    if((f(a) * f(b)) > 0) {
        cout<<"No hay raiz o hay un numero par de ellas"<<endl;
        exit(0);
    }
    double error = 1;
    int iteration = 0;

    do {
        root = ((a+b)/2);
        iteration++;
    if((f(a) * f(root)) > 0) {
        a = root;
    } else if((f(a) * f(root)) < 0) {
        b = root;
    } else {
        cout<<"La raiz es "<<root<<endl;
        exit(0);
    }
    error = ((b-a)/2);
    
    } while(error > tolerance);

    //Variamos el %.f dependiendo cuanto nos de el error, en base a ese error redondeamos a la misma cifra significativa a la raiz
    printf("La raiz es %.5f\n", root);
    printf("El error es %.5f\n", error);
    printf("Raiz: %.5f +- %.5f\n", root, error);
    printf("Se ha iterado %d vez/veces\n", iteration);

    return 0;
}

double f(double x) {
    return (-2+(7*x)-(5*pow(x, 2))+(6*pow(x, 3)));
}

