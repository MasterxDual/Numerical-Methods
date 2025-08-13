/*La eficiencia de este motor para gases monoatómi os está dada por una ecuacion donde pasé el 0.3 restando al otro lado y luego reemplace los valores
en la ecuación para que me quede como la función del programa.

γ = 5/3 . Encontrar T2/T1 tal que la eficiencia sea del 30 % (η = 0,3).

a = 5.3, b = 5.5, tolerance = 1e-3;

Este resultado me dió con el método común (Exercise3Biseccion.cpp)
La raiz es 5.414
El error es 0.001
Raiz: 5.414 +- 0.001
Se ha iterado 7 vez/veces

Este resultado me dió con el método modificado que está en este programa
La raiz es 5.413
El error es 0.001
Raiz: 5.413 +- 0.001
Se ha iterado 8 vez/veces

T2/T1 es igual a 5.414 +-0.001
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
    double old_root = a;

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
    error = fabs(root-old_root);
    old_root = root;
    
    } while(error > tolerance);
    
    //Variamos el %.f dependiendo cuanto nos de el error, en base a ese error redondeamos a la misma cifra significativa a la raiz
    printf("La raiz es %.3f\n", root);
    printf("El error es %.3f\n", error);
    printf("Raiz: %.3f +- %.3f\n", root, error);
    printf("Se ha iterado %d vez/veces\n", iteration);

    return 0;
}

double f(double x) {
    return (((log(x)-(1-pow(x, -1)))/(log(x)+((3*(1-pow(x, -1)))/2)))-0.3);
}