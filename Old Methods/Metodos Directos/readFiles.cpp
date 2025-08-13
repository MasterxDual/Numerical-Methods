#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <iostream>
using namespace std;


int main(int argc, char* argv[])
{
    FILE *fp;
    char ch;
    fp = fopen("datos.dat","r");
    if ( fp == NULL )
    {
        puts ( "No se puede abrir el archivo");
    }

    //contador de filas
    int filas=0;
    char c;
    int maxValues = 0;
	int columnas;
    
    while((c = fgetc(fp)) != EOF)
    {
        if(c == '\n'){
            filas++;
		}
    }

    printf("numero de filas=%i\n",filas);

    //cerramos y reabrimos el archivo para resetear el puntero
    fclose(fp);
    fp = fopen("datos.dat","r");
    double m[filas][100];


    //Cargo los datos leidos en el array
    int i, j;
	for(i = 0; i < filas; i++) {
		j = 0;
		do {
        fscanf(fp, "%lf", &(m[i][j]));
        j++;
			} while((c = fgetc(fp)) != '\n');
	}

	columnas = j;
	
	printf("numero de columnas= %i\n",columnas);
    //imprimo la matriz para verificar que lo leyo correctamente
    
        for(i=0;i<filas;i++)
    {
        for(j=0;j<columnas;j++)
        {
            printf("%lf ",m[i][j]);
        }
        printf("\n");
    }
}

