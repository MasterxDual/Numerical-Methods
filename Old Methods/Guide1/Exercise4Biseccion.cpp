/*Encontrar el número n de bisecciones sucesivas que garanticen que el punto medio Cn es una
aproximación a un cero (de la ecuación considerada) con un error menor que un valor prefijado
δ. Tenga en cuenta que dado dos extremos de intervalo a y b, encontrado el punto medio c, y
considerando un cero r de la función, se cumple en una bisección la siguiente desigualdad:
|r − c| ≤ |b − a|/ 2

Partiendo de esta desigualdad, encuentre el valor de n para una tolerancia δ = |r − c| requerida.

Sabemos que, cada vez que hacemos mas iteraciones, el error se divide en dos, quedandonos resumidamente
que error = b-a/2**n (2 elevado a la n), ya que dividimos en dos la ecuación mencionada sucesivamente a medida que iteramos. 
Por lo tanto, con esta ecuación, nos queda que b-a/2**n = δ y, despejando n de esta ecuación, nos queda que 
n = ln(b-a/δ) / ln(2).

El resto del programa esta de adorno, no hay que utilizar ningun metodo matematico.
*/

#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <stdio.h>

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