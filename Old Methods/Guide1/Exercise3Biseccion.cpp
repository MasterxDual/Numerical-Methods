/*Calcule la masa m necesaria del paracaidista tal que la velocidad de caída sea v = 35 m
seg en t = 7 seg. Utilice el método de bisección y de regula falsi.
 
a = 50, b = 70, tolerance = 1e-3
La raiz es 63.649
El error es 0.001
Raiz: 63.649 +- 0.001
Se ha iterado 14 vez/veces

La masa m necesaria es igual a 63.649 kilogramos.
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
    printf("La raiz es %.3f\n", root);
    printf("El error es %.3f\n", error);
    printf("Raiz: %.3f +- %.3f\n", root, error);
    printf("Se ha iterado %d vez/veces\n", iteration);

    return 0;
}

double f(double x) {
    return ((0.7*x)*(1-exp(-98/x))-35);
}