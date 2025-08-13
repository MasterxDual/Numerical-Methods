#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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

    /*Lo unico que varía en este método (con respecto al primer método) 
    es la fórmula para calcular la raiz, aquí se usa una recta tangente a los dos puntos que intersectan a la funcion en el eje y.*/

    do {
        root = ((a * f(b)) - (b * f(a))) / (f(b) - f(a));
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
    return (log(x)+exp(sin(x))-x);
}