/*Consideremos el problema de hallar la porción de una esfera de radio r que queda sumergida
(parcialmente) en agua. Supongamos que la esfera está construida en madera de pino con una
densidad ρ = 0,638 gr/cm3 y que su radio mide r = 10 [cm]. Cual es la profundidad d a la que
está sumergida la esfera?. El volumen de la esfera es Ve = 4/3 * π *r3 .

Sabemos que ρ = m/V, y Ve = 4/3 * π *r3, por lo tanto Ve = 4/3 * π * (10)al cubo osea tenemos que Ve = 4188.79020479 cm cubicos.
Luego despejamos la masa de la formula de densidad, quedando que m = ρ * Ve, reemplazando en esta ecuacion nos queda que 
m = 2672.44802 gramos. Luego por Fuerzas hidrostaticas la esfera esta flotando (en equilibrio) en el agua, por lo tanto tenemos que
Pe - PH2O = 0, ademas PH20 = ρH2O * g * VpedazoEsferaSumergido y Pe = me * g. Ademas, sabemos que la altura de la esfera que no esta dentro 
del agua mas la altura de la esfera que esta dentro del agua tiene que ser igual al doble del radio, es decir, h + d = 2*r (esta ecuacion no fue utilizada). 
Nos faltaria saber a que es igual el VpedazoEsferaSumergido. Para esto buscamos en Internet y sabemos que el volumen de un casquete esferico es igual a:
VpedazoEsferaSumergido = VcasqueteEsferico = ((pi * d) / 6) * (3*a**2 + h**2) o tambien VcasqueteEsferico = (pi * d**2 / 3) * (3*r - d) y utilizaremos esta ultima formula.
Sabiendo que Pe = PH2O, reemplazando en esta igualdad nos queda que me * g = ρH2O * g * VpedazoEsferaSumergido, reemplazando el
VpedazoEsferaSumergido = (pi * d**2 / 3) * (3*r - d) en la igualdad nos queda que me * g = ρH2O * g * ((pi * d**2 / 3) * (3*r - d)) y, reordenando en esta ecuacion
para despejar d que es la incognita que nos pide el ejercicio, nos queda una cuadratica igualada a 0 -> ...
... -> ρH2O * g * pi * r * d**2 - ρH2O * g * (pi / 3) * d - (me * g)) = 0 donde viendono en terminos de x nuestra incognita y reemplazando las constantes
nos queda -> 30809.6 (gr/(cm * s**2)) * x**2 - 1026.99 (gr/(cm**2 * s**2)) * x - 2620861.9 ((gr * cm) / s**2) = 0.
Utilicé el método de Newton-Raphson ya que la derivada es sencilla de conseguir. Por intuicion se probó que es mas sencillo utilizar este método antes que 
el de punto fijo, ya que para el de punto fijo debería despejar la x, y haciendo eso nos quedaría una función mas compleja de utilizar. En cambio con el de 
Newton-Raphson con solo derivar fue suficiente.

X0 = 8, tolerance = 1e-3;
La primer raiz de f(x) es igual a 9.23992
El error es igual a 0.00050
Resultado : 9.23992 +- 0.00050
Se han iterado 3 veces

X0 = -8, tolerance = 1e-3;
La segunda raiz de f(x) es igual a -9.20658
El error es igual a 0.00044
Resultado : -9.20658 +- 0.00044
Se han iterado 3 veces

Viendo estos resultados, concluimos que la raiz que tiene sentido físico es la positiva (primer raiz), por lo que la distancia o profundidad d a la que está 
sumergida la esfera es igual a 9.24 centimetros -> d = 9.24 cm.
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

    printf("La raiz de f(x) es igual a %.5f\n", newX);
    printf("El error es igual a %.5f\n", error); 
    printf("Resultado : %.5f +- %.5f\n", newX, error);
    printf("Se han iterado %d veces\n", iteration);

    return 0;
}

double fPrima(double x) {
    return (61619.2*x - 1026.99);
}

double f(double x) {
    return (30809*pow(x, 2) - 1026.99*x - 2620861.9);
}


