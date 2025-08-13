/*Este ejercicio se realizó tomando una forma rectangular.

Se construye un contenedor sin tapa a partir de una hoja metálica rectangular que mide 10 ×
16 [cm].Cuál debe ser el lado de los cuadrados que hay que recortar en cada esquina para que el
volumen del contenedor sea 100 cm3 ?. Precisión 1 × 10−9 [cm]. Note que puede formar más de
una geometría de contenedor con este requerimiento, alguna de estas pueden no tener sentido físico,
verifique. 

Utilicé el metodo de Newton-Raphson ya que viendo como quedaría la derivada de f(x) evaluada en algun 
punto Xi, no nos da ni muy grande ni muy chico el resultado, por lo que este metodo es perfecto.

Sabemos que el area del rectangulo debe ser igual a 16 * 10 = 160 cm2, y, viendo una grafica hecha 
a mano, tenemos que el area del contenedor debe ser igual a (16 - (2*x)) * (10 - (2*x)) = Ac, y 
de aqui podemos calcular el volumen del contenedor que es igual a Vc = ((16 - (2*x)) * (10 - (2*x))) * x
donde x es igual al lado de cualquiera de los 4 cuadrados que hay que recortar en las esquinas.
Resolviendo la ecuacion de Vc nos queda que Vc = (4 * xcubo) - (52 * xcuadrado) + (160 * x) y esto debe
ser igual a 100 cmcubicos osea 100 = (4 * xcubo) - (52 * xcuadrado) + (160 * x) y dejando esta 
ecuacion igualada a 0 nos queda que (4 * xcubo) - (52 * xcuadrado) + (160 * x) - 100 = 0. De aqui utilizamos
el metodo de Newton-Raphson y nos da las siguientes raices, que serian los lados de los cuadrados 
que nos pide el problema.

X0 = 0, tolerance = 1e-9;
La primer raiz de f(x) es igual a 0.8390188831
El error es igual a 0.0000000000
Resultado : 0.8390188831 +- 0.0000000000
Se han iterado 6 veces

X0 = 3, tolerance = 1e-9;
La segunda raiz de f(x) es igual a 3.4017486475
El error es igual a 0.0000000000
Resultado : 3.4017486475 +- 0.0000000000
Se han iterado 5 veces

X0 = 8, tolerance = 1e-9;
La tercer raiz de f(x) es igual a 8.7592324694
El error es igual a 0.0000000000
Resultado : 8.7592324694 +- 0.0000000000
Se han iterado 6 veces

El lado de los cuadrados debe ser igual a 0.839 cm, o puede ser 3.401 cm, o puede ser 8.759 cm.
*/

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

    printf("La raiz de f(x) es igual a %.10f\n", newX);
    printf("El error es igual a %.10f\n", error); 
    printf("Resultado : %.10f +- %.10f\n", newX, error);
    printf("Se han iterado %d veces\n", iteration);

    return 0;
}

double fPrima(double x) {
    return ((12 * pow(x, 2)) - (104 * x) + 160);
}

double f(double x) {
    return ((4 * pow(x, 3)) - (52 * pow(x, 2)) + (160 * x) - 100);
}
