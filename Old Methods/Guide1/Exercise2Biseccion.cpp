/*Obtenga las raíces de la siguiente función utilizando los métodos de bisección y de régula falsi:
g (a) = a10 − 1 

1. Verifique la convergencia de ambos métodos para un error ε = 1 × 10−5. Cuantas iteraciones
le tomó a cada uno de los métodos encontrar las raíces con la precisión deseada?. Obtenga
una gráfica comparativa del error de aproximación de cada método en función del número de
iteraciones.

Metodo de Bisección
a = -1.05, b = -0.9, tolerance = 1e-5;
La primer raiz es -1.000000
El error es 0.000007
Primer Raiz: -1.000000 +- 0.000007
Se ha iterado 21 vez/veces

a = 0.9, b = 1.2, tolerance = 1e-5;
La segunda raiz es 1.000000
El error es 0.000007
Segunda Raiz: 1.000000 +- 0.000007
Se ha iterado 22 vez/veces.

Comparando ambos métodos, al método de Regula Falsi le tomó menos iteraciones que al de bisección, esto demuestra que es mas efectivo (solo en los casos que
las funciones analizadas no crezcan rapidamente cerca de sus raices).
En general, el método de la regla falsa tiende a converger más rápido que el método de bisección.
Esto se debe a que la regla falsa utiliza una interpolación lineal para estimar la raíz en cada iteración, 
lo que puede proporcionar aproximaciones más precisas y permitir una convergencia más rápida en funciones que cambian rápidamente.
*/

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;

double g(double a);

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

    if((g(a) * g(b)) > 0) {
        cout<<"No hay raiz o hay un numero par de ellas"<<endl;
        exit(0);
    }
    double error = 1;
    int iteration = 0;
    double old_root = a; //Nueva variable para comparar los errores de aproximacion

    //Escribir un archivo con los datos
    ofstream ofs; //Creo un objeto de la libreria ofstream
    ofs.open("TablaBiseccion.txt"); //Abro un archivo, si no existe se crea
    ofs<<"Metodo Biseccion\nIteracion\tError\n"<<endl; //Escribo en el archivo

    do {
        root = ((a+b)/2);
        iteration++;
    if((g(a) * g(root)) > 0) {
        a = root;
    } else if((g(a) * g(root)) < 0) {
        b = root;
    } else {
        cout<<"La raiz es "<<root<<endl;
        exit(0);
    }
    //Ahora utilizo este método para graficar las funciones comparando los errores de aproximacion
    error = fabs((root - old_root) / root) * 100;
    ofs<<fixed<<setprecision(6)<<iteration<<"\t\t"<<fixed<<setprecision(6)<<error<<endl; //Escribo en el archivo creado
    old_root = root;
    } while(error > tolerance);
    ofs.close(); 

    cout<<"LECTURA DE LA TABLA CREADA"<<endl;
    string message = " ";
    ifstream ifs; //Creo un objeto de la libreria ifstream
    ifs.open("TablaBiseccion.txt");
    while(!ifs.eof()) { //Recorre el archivo hasta que llega al fin de linea
        getline(ifs, message);
        cout<<message<<endl;
    }
    ifs.close();

    //Variamos el %.f dependiendo cuanto nos de el error, en base a ese error redondeamos a la misma cifra significativa a la raiz
    printf("La primer raiz es %.6f\n", root);
    printf("El error es %.6f\n", error);
    printf("Primer Raiz: %.6f +- %.6f\n", root, error);
    printf("Se ha iterado %d vez/veces\n", iteration);
    return 0;
}

double g(double a) {
    return (pow(a, 10) - 1);
}