/*Sea g (x) = x2 + x − 4. Podemos utilizar iteración de punto fijo para hallar las soluciones de la
ecuación x = g (x)?. Por que?. Suponiendo una ecuación arbitraria que tiene un punto fijo P , Cual
es la ventaja de tener g'(P ) ≈ 0 en un proceso de iteración de punto fijo?.

Primer intento:
Haciendo la derivada de g(x) tenemos que g'(x) = 2*x + 1, luego fabs(g'(x)) > 1 por lo que nunca converge hacia el punto fijo real de g(x). 
Es decir, en cada iteración me alejo más del punto fijo real, haciendo que el programa no encuentre nunca el punto fijo (valor de x donde corta la recta
y = x con la funcion g(x), siendo este punto fijo la raiz de f(x)).

Depende como sea hayamos tomado a la funcion dada, en este caso g(x) lo tomamos que es f(x) y luego conseguimos g(x) para hacer luego su derivada -> {
Segundo intento:
Despejamos x:
Tenemos que f(x) = x2 + x -4 = 0, luego x2 - 4 = -x, o equivalentemente x = 4 - x2, viendo que g(x) = 4 - x2. Entonces derivamos a g(x), teniendo
g'(x) = -2*x, luego fabs(g'(x)) = 2 * fabs(x) teniendo que fabs(g'(x)) > 1, llegando a la conclusion de que tampoco converge en ningun momento.

Tercer intento: 
Sumamos x a ambos miembros:
Nos queda x2 + 2*x - 4 = x, luego g(x) = x2 + 2*x -4, luego g'(x) = 2*x + 2, es decir g'(x) = 2*(x + 1), aplicando el valor absoluto nos queda
fabs(g'(x)) = 2 * fabs(x + 1), y esto quiere decir que fabs(g'(x)) > 1. Concluyendo que tampoco converge en ningun momento.
}

Si suponemos una ecuacion arbritaria que tiene un punto fijo P, la ventaja de tener g'(P) ≈ 0 quiere decir que g'(P) < 1, lo que nos quiere decir que la derivada 
evaluada en el punto fijo P nos dará un valor menor a 1, lo que nos dice que en cada iteracion nos hemos acercado cada vez mas al punto fijo P y ademas
hemos podido encontrar el punto fijo y que el programa pudo realizar
todos los calculos necesarios para llegar a la raiz de f(x) o al punto fijo de g(x), habiendo convergido la función. 
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

    printf("El punto fijo de g(x) o la raiz de f(x) es igual a %.5f\n", newX);
    printf("El error es igual a %.5f\n", error); 
    printf("Resultado : %.5f +- %.5f\n", newX, error);
    printf("Se han iterado %d veces\n", iteration);

    return 0;
}

double g(double x) {
    //return (pow(x, 2) + x - 4); Funcion del primer intento que esta detallado arriba
    return (4 - pow(x, 2)); // Funcion del segundo intento que esta detallado arriba
}

double f(double x) {
    return ((1/2) * sin(x) + (2*x));
}