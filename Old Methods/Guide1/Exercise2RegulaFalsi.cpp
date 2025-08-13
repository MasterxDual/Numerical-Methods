/*Obtenga las raíces de la siguiente función utilizando los métodos de bisección y de régula falsi:
g (a) = a10 − 1 

1. Verifique la convergencia de ambos métodos para un error ε = 1 × 10−5. Cuantas iteraciones
le tomó a cada uno de los métodos encontrar las raíces con la precisión deseada?. Obtenga
una gráfica comparativa del error de aproximación de cada método en función del número de
iteraciones.

Metodo Regula Falsi
a = -1.05, b = -0.9, tolerance = 1e-5;
La primer raiz es -1.00
El error es 0.00
Primer Raiz: -1.00 +- 0.00
Se ha iterado 10 vez/veces

a = 0.9, b = 1.2, tolerance = 1e-5;
La segunda raiz es 1.00
El error es 0.00
Segunda Raiz: 1.00 +- 0.00
Se ha iterado 28 vez/veces

En la segunda raiz, el metodo de Regula Falsi contiene una mayor cantidad de iteraciones que el metodo de Biseccion ya que en esa raiz, la funcion crece
rapidamente, haciendo que disminuya la convergencia del metodo de Regula Falsi.
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

    double error = 1;
    int iteration = 0;
    double old_root = a; //Nueva variable para comparar los errores de aproximacion

    //Escribir un archivo con los datos
    ofstream ofs; //Creo un objeto de la libreria ofstream
    ofs.open("TablaRegulaFalsi.txt"); //Abro un archivo, si no existe se crea
    ofs<<"Metodo Regula Falsi\nIteracion\tError\n"<<endl; //Escribo en el archivo

    if(g(a) * g(b) < 0) {
    do {
        root = ((a * g(b)) - (b * g(a))) / (g(b) - g(a));
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
    ifs.open("TablaRegulaFalsi.txt");
    while(!ifs.eof()) { //Recorre el archivo hasta que llega al fin de linea
        getline(ifs, message);
        cout<<message<<endl;
    }
    ifs.close();

    //Variamos el %.f dependiendo cuanto nos de el error, en base a ese error redondeamos a la misma cigra signigicativa a la raiz
    printf("La primer raiz es %.2f\n", root);
    printf("El error es %.2f\n", error);
    printf("Primer Raiz: %.2f +- %.2f\n", root, error);
    printf("Se ha iterado %d vez/veces\n", iteration);
    } else {
        cout<<"No hay raiz o hay un numero par de ellas"<<endl;
    }
    return 0;
}

double g(double a) {
    return (pow(a, 10) - 1);
}