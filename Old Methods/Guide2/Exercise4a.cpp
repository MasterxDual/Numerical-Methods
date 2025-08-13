/*Este ejercicio se realizó tomando una forma cilindrica, que no es la que debia ser. De todos modos
la dejo aquí plasmada para mostrar como se haría
Se construye un contenedor sin tapa a partir de una hoja metálica rectangular que mide 10 ×
16 [cm].Cuál debe ser el lado de los cuadrados que hay que recortar en cada esquina para que el
volumen del contenedor sea 100 cm3 ?. Precisión 1 × 10−9 [cm]. Note que puede formar más de
una geometría de contenedor con este requerimiento, alguna de estas pueden no tener sentido físico,
verifique.


Utilicé el metodo de Newton-Raphson ya que viendo como quedaría la derivada de f(x) evaluada en algun 
punto Xi, no nos da ni muy grande ni muy chico el resultado, por lo que este metodo es perfecto.

Con un X0 = 1 y una tolerancia = 1e-9 nos da:
La primer raiz de f(x) es igual a 1.29238
El error es igual a 0.00000
Resultado : 1.29238 +- 0.00000
Se han iterado 5 veces
Es decir, r1 = 1.29238 cm.

Con un X0 = 7 y una tolerance = 1e-9 nos da:
La segunda raiz de f(x) es igual a 6.40199
El error es igual a 0.00000
Resultado : 6.40199 +- 0.00000
Se han iterado 5 veces
Es decir, r2 = 6.40199 cm.

Con un X0 = -8 y una tolerance = 1e-9 nos da:
La tercer raiz de f(x) es igual a -7.69438
El error es igual a 0.00000
Resultado : -7.69438 +- 0.00000
Se han iterado 3 veces
Es decir, r3 = -7.69438 cm.

Pero como la ecuacion dice que el Area del cilindro debe ser MENOR o igual, podemos bajar los valores de r1, r2 y r3 (ya que probé con esos valores exactos y no 
me da el calculo final). Por lo tanto usaremos r1 = 0.9 cm, r2 = 5 cm y r3 = -9 cm.

Reemplazando estas raices en h = 100/pi*r2 nos queda:
h1 = 25.7831007809 cm, h2 = 795.774715459cm y h3 = 2578.31007809 cm

Reemplazando estas alturas y estos radios en Ac = pi*r2 + 2*pi*r*h
nos queda que Ac1 = 148.344690049 cm2 , Ac2 = numero muy grande y Ac3 = numero mas grande todavia (estos dos ultimos no cumplen con lo que pide el problema)
Viendo estos resultados, deberia seguir aproximando h2 y h3 a valores menores, pero con r1 no pasa eso, lo aproxime a un numero bastante cercano a su valor original
por lo que fisicamente posible es con r1 solamente.

Finalmente tenemos que Ac1 <= 160 cm2. Pero esto no terminó todavía. 

El area de un cuadrado es igual a lado al cuadrado, osea Acuadrado = l2, necesitamos saber el lado de alguno de los 4 cuadrados que recortaremos en cada
esquina del cilindro, entonces tenemos que el Areatransversal es igual a 2*pi*r*h, osea Atransv = 145.782465497 cm2, y de aquí debemos realizar el ultimo calculo:
Primero con Pitágoras calculamos la diagonal de un cuadrado, es decir, la hipotenusa, que es igual a H = raizcuadrada(lado2 + lado2) = raizcuadrada(2)*lado
Luego el perimetro de la base del cilindro es igual a P = 2*pi*r, luego divido al perimetro entre cuatro ya que hay cuatro cuadrados, uno en
cada esquina, P = pi*r/2, luego igualamos la diagonal del cuadrado a la fraccion del perimetro:
P = H -> pi*r/2 = raizcuadrada(2)*lado y luego despejamos el lado que nos queda -> lado = pi*r/raizcuadrada(2)*2. 
Reemplazando r1 = 0.9 nos queda que lado =  0.99964866108 cm. Ver esto ultimo de los cuadrados
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
    return (2*M_PI*x - (200/pow(x, 2)));
}

double f(double x) {
    return (M_PI*pow(x, 2) + (200/x) - 160);
}

/*A f(x) la conseguí de la siguiente manera: 
Sabemos que el el area de la hoja metalica es igual a: Areah = 10*16 = 160 cm2. Ademas, el volumen del contenedor sin tapa (cilindro sin el circulo de arriba)
es igual a: Vc = pi*r2*h = Areabase * h donde h y r son incognitas. Luego tenemos, además, que el volumen del contenedor debe ser igual a 100 cm3. Es decir, tenemos que
Vc = 100 cm3. Seguimos con formulas, El area del contenedor debe ser igual a: Ac = pi*r2 + 2*pi*r*h = Arealateral + Areabase ya que tenemos un solo circulo mas el area 
lateral del cilindro. Luego tenemos que el Area del contenedor debe ser menor o igual al Area de la hoja metalica, que eso es lo que nos dice el ejercicio 
implícitamente. Nos queda que Ac <= Areah es decir pi*r2 + 2*pi*r*h <= 160, luego de la ecuacion de Vc = 100 cm3, reemplazamos y nos queda pi*r2*h = 100, despejamos 
h y nos queda que h = 100/(pi*r2), reemplazando esta ecuacion en pi*r2 + 2*pi*r*h <= 160, y luego dejando el cero del lado derecho de esta ecuacion reemplazada nos 
queda pi*r2 + 200/r - 160 <= 0, aproximando a cero esta ecuacion para encontrar las raices con algun metodo numerico nos queda pi*r2 + 200/r - 160 = 0, encontrando 
las raices de esta ecuacion y asi conseguimos a f(x) y luego derivando
a f(x) nos quedó fPrima(x) que estan arriba de este texto.*/