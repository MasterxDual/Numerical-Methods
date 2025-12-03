/**
 * ============================================================================
 * PROBLEMA N° 1: Interpolación con Splines e Integración
 * ============================================================================
 * 
 * Función: f(x) = e^(2·sin(x)) · [1 + ln(x)]
 * 
 * TAREAS:
 * a) Generar tabla equiespaciada con h=0.1 usando splines cúbicas
 *    - Calcular valores exactos y error
 * b) Calcular ∫[1,2] f(x)dx con Simpson 1/3 usando valores del spline
 * c) Calcular ∫[1,2] f(x)dx con Simpson 1/3 usando función exacta
 * 
 * ============================================================================
 Para usarlo:
    * 1. Asegúrate de tener un archivo "nodos.txt" con los puntos de datos.
    * 2. Compila el programa usando un compilador compatible con C17.
    * 3. Ejecuta el programa y sigue las instrucciones en pantalla.
 * ============================================================================
 
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* ============================================================================
   CONFIGURACIÓN DE LA FUNCIÓN A ANALIZAR
   ============================================================================
   IMPORTANTE: Cuando cambies la función f_exacta(), actualiza también esta
   cadena para que coincida con la expresión matemática.
   ============================================================================ */

#define FUNCION_DESCRIPCION "f(x) = e^(2·sin(x)) · [1 + ln(x)]"

/* ============================================================================
   PROTOTIPOS DE FUNCIONES
   ============================================================================ */

double f_exacta(double x);
void construirSistemaSplines(double *x, double *y, int n, double *A, double *b_vec);
void eliminacionGaussiana(double *A, double *b_vec, double *solution, int n);
double evaluarSpline(double x_eval, double *x, double *y, double *M, int n);
void generarTablaEquiespaciada(double *x_orig, double *y_orig, int n_orig, 
                               double **x_equi, double **y_equi, int *n_equi, double h_objetivo);
int leerDatosDesdeArchivo(const char *filename, double **x, double **y);

/* ============================================================================
   FUNCIÓN PRINCIPAL
   ============================================================================ */

int main(void)
{
    /* Leer datos originales desde el archivo nodos.txt */
    double *x_original = NULL;
    double *y_original = NULL;
    int n_original = 0;
    
    printf("\n╔═════════════════════════╗\n");
    printf("║   Splines Cúbicas       ║\n");
    printf("╚═════════════════════════╝\n");
    
    printf("\nLeyendo datos desde 'nodos.txt'...\n");
    n_original = leerDatosDesdeArchivo("nodos.txt", &x_original, &y_original);
    
    if (n_original == 0) {
        printf("[ERROR] No se pudieron leer los datos del archivo.\n");
        return 1;
    }
    
    printf("✓ Se leyeron %d puntos correctamente\n", n_original);
    printf("\nDatos originales:\n");
    printf("  Intervalo: [%.2lf, %.2lf]\n", x_original[0], x_original[n_original-1]);
    printf("  Primer punto: (%.2lf, %.3lf)\n", x_original[0], y_original[0]);
    printf("  Último punto: (%.2lf, %.3lf)\n", x_original[n_original-1], y_original[n_original-1]);
    
    double *x_equi = NULL;
    double *y_equi = NULL;
    int n_equi = 0;
    
    printf("\nFunción: %s\n", FUNCION_DESCRIPCION);
    
    /* ========================================================================
       SOLICITAR EL ESPACIADO h AL USUARIO
       ======================================================================== */
    
    double h_objetivo;
    printf("\n╔════════════════════════════════════════════════════════════════╗\n");
    printf("║  CONFIGURACIÓN DE PARÁMETROS                                   ║\n");
    printf("╚════════════════════════════════════════════════════════════════╝\n");
    printf("\nIngrese el espaciado h deseado para la tabla equiespaciada: ");
    scanf("%lf", &h_objetivo);
    while (getchar() != '\n');  // Limpiar buffer
    
    if (h_objetivo <= 0) {
        printf("[ERROR] El espaciado h debe ser positivo.\n");
        free(x_original);
        free(y_original);
        return 1;
    }
    
    double rango = x_original[n_original-1] - x_original[0];
    int puntos_estimados = (int)(rango / h_objetivo) + 1;
    printf("  ✓ Con h = %.4lf se generarán aproximadamente %d puntos\n", h_objetivo, puntos_estimados);
    
    /* ========================================================================
       GENERAR TABLA EQUIESPACIADA CON h ESPECIFICADO
       ======================================================================== */
    
    printf("\n╔══════════════════════════════════════════════════════╗\n");
    printf("║   Tabla Equiespaciada con Splines (h = %.4lf)       ║\n", h_objetivo);
    printf("╚══════════════════════════════════════════════════════╝\n");
    generarTablaEquiespaciada(x_original, y_original, n_original,
                             &x_equi, &y_equi, &n_equi, h_objetivo);
    
    /* Calcular valores exactos y errores */
    printf("\n┌──────┬──────────┬────────────────┬────────────────┬──────────────┐\n");
    printf("│  i   │   x_i    │  y (spline)    │  y (exacta)    │   Error      │\n");
    printf("├──────┼──────────┼────────────────┼────────────────┼──────────────┤\n");
    
    double error_max = 0.0;
    double error_prom = 0.0;
    
    for (int i = 0; i < n_equi; i++) {
        double y_exacta = f_exacta(x_equi[i]);
        double error = fabs(y_equi[i] - y_exacta);
        
        error_prom += error;
        if (error > error_max) error_max = error;
        
        printf("│ %3d  │  %.2lf    │  %12.6lf  │  %12.6lf  │  %10.6lf  │\n",
               i, x_equi[i], y_equi[i], y_exacta, error);
    }
    
    printf("└──────┴──────────┴────────────────┴────────────────┴──────────────┘\n");
    
    error_prom /= n_equi;
    
    printf("\nESTADÍSTICAS DE ERROR:\n");
    printf("  Error máximo:   %.6lf\n", error_max);
    printf("  Error promedio: %.6lf\n", error_prom);
    

    
    printf("\nCONCLUSIONES:\n");
    printf("  • Las splines cúbicas aproximan bien la función\n");

    /* ========================================================================
       GUARDAR RESULTADOS
       ======================================================================== */
    
    FILE *archivo = fopen("resultados.txt", "w");
    if (archivo != NULL) {
        fprintf(archivo, "# Splines\n");
        fprintf(archivo, "# %s\n\n", FUNCION_DESCRIPCION);
        
        fprintf(archivo, "# Tabla equiespaciada (h = %.4lf)\n", h_objetivo);
        fprintf(archivo, "# i\tx\ty_spline\ty_exacta\terror\n");
        for (int i = 0; i < n_equi; i++) {
            fprintf(archivo, "%d\t%.2lf\t%.6lf\t%.6lf\t%.6lf\n",
                    i, x_equi[i], y_equi[i], f_exacta(x_equi[i]), 
                    fabs(y_equi[i] - f_exacta(x_equi[i])));
        }

        fclose(archivo);
        printf("\n✓ Resultados guardados en 'resultados.txt'\n");
    }
    
    /* Liberar memoria */
    free(x_original);
    free(y_original);
    free(x_equi);
    free(y_equi);
    
    return 0;
}

/* ============================================================================
   IMPLEMENTACIÓN DE FUNCIONES
   ============================================================================ */

/**
 * Función exacta: f(x) = e^(2·sin(x)) · [1 + ln(x)]
 */
double f_exacta(double x)
{
    return exp(2.0 * sin(x)) * (1.0 + log(x));
}

/**
 * Genera tabla equiespaciada usando splines cúbicas
 */
void generarTablaEquiespaciada(double *x_orig, double *y_orig, int n_orig,
                               double **x_equi, double **y_equi, int *n_equi, double h_objetivo)
{
    printf("\nConstruyendo splines cúbicas...\n");
    
    /* Construir sistema para splines */
    double *A = (double *)malloc(n_orig * n_orig * sizeof(double));
    double *b_vec = (double *)malloc(n_orig * sizeof(double));
    double *M = (double *)malloc(n_orig * sizeof(double));
    
    construirSistemaSplines(x_orig, y_orig, n_orig, A, b_vec);
    eliminacionGaussiana(A, b_vec, M, n_orig);
    
    printf("Splines calculadas\n");
    
    /* Generar tabla equiespaciada */
    double a = x_orig[0];
    double b = x_orig[n_orig - 1];
    *n_equi = (int)((b - a) / h_objetivo) + 1;
    
    *x_equi = (double *)malloc(*n_equi * sizeof(double));
    *y_equi = (double *)malloc(*n_equi * sizeof(double));
    
    printf("Generando %d puntos equiespaciados...\n", *n_equi);
    
    for (int i = 0; i < *n_equi; i++) {
        (*x_equi)[i] = a + i * h_objetivo;
        (*y_equi)[i] = evaluarSpline((*x_equi)[i], x_orig, y_orig, M, n_orig);
    }
    
    printf("Tabla generada\n");
    
    free(A);
    free(b_vec);
    free(M);
}

/**
 * Construye sistema tridiagonal para splines cúbicas naturales
 */
void construirSistemaSplines(double *x, double *y, int n, double *A, double *b_vec)
{
    for (int i = 0; i < n * n; i++) A[i] = 0.0;
    
    double *h = (double *)malloc((n-1) * sizeof(double));
    for (int i = 0; i < n-1; i++) {
        h[i] = x[i+1] - x[i];
    }
    
    A[0] = 1.0;
    b_vec[0] = 0.0;
    
    for (int i = 1; i < n-1; i++) {
        A[i * n + (i-1)] = h[i-1];
        A[i * n + i] = 2.0 * (h[i-1] + h[i]);
        A[i * n + (i+1)] = h[i];
        b_vec[i] = 6.0 * ((y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1]);
    }
    
    A[(n-1) * n + (n-1)] = 1.0;
    b_vec[n-1] = 0.0;
    
    free(h);
}

/**
 * Eliminación Gaussiana con pivoteo parcial
 */
void eliminacionGaussiana(double *A, double *b_vec, double *solution, int n)
{
    double *A_copy = (double *)malloc(n * n * sizeof(double));
    double *b_copy = (double *)malloc(n * sizeof(double));
    
    for (int i = 0; i < n * n; i++) A_copy[i] = A[i];
    for (int i = 0; i < n; i++) b_copy[i] = b_vec[i];
    
    for (int k = 0; k < n-1; k++) {
        int max_row = k;
        double max_val = fabs(A_copy[k * n + k]);
        
        for (int i = k+1; i < n; i++) {
            if (fabs(A_copy[i * n + k]) > max_val) {
                max_val = fabs(A_copy[i * n + k]);
                max_row = i;
            }
        }
        
        if (max_row != k) {
            for (int j = 0; j < n; j++) {
                double temp = A_copy[k * n + j];
                A_copy[k * n + j] = A_copy[max_row * n + j];
                A_copy[max_row * n + j] = temp;
            }
            double temp = b_copy[k];
            b_copy[k] = b_copy[max_row];
            b_copy[max_row] = temp;
        }
        
        for (int i = k+1; i < n; i++) {
            double factor = A_copy[i * n + k] / A_copy[k * n + k];
            for (int j = k; j < n; j++) {
                A_copy[i * n + j] -= factor * A_copy[k * n + j];
            }
            b_copy[i] -= factor * b_copy[k];
        }
    }
    
    for (int i = n-1; i >= 0; i--) {
        solution[i] = b_copy[i];
        for (int j = i+1; j < n; j++) {
            solution[i] -= A_copy[i * n + j] * solution[j];
        }
        solution[i] /= A_copy[i * n + i];
    }
    
    free(A_copy);
    free(b_copy);
}

/**
 * Evalúa spline cúbica en un punto
 */
double evaluarSpline(double x_eval, double *x, double *y, double *M, int n)
{
    int j = 0;
    for (int i = 0; i < n-1; i++) {
        if (x_eval >= x[i] && x_eval <= x[i+1]) {
            j = i;
            break;
        }
    }
    
    if (x_eval < x[0]) j = 0;
    if (x_eval > x[n-1]) j = n-2;
    
    double h_j = x[j+1] - x[j];
    
    double t1 = (M[j] / 6.0) * pow(x[j+1] - x_eval, 3) / h_j;
    double t2 = (M[j+1] / 6.0) * pow(x_eval - x[j], 3) / h_j;
    double t3 = (y[j] - M[j] * h_j * h_j / 6.0) * (x[j+1] - x_eval) / h_j;
    double t4 = (y[j+1] - M[j+1] * h_j * h_j / 6.0) * (x_eval - x[j]) / h_j;
    
    return t1 + t2 + t3 + t4;
}

/**
 * Lee datos desde archivo de texto
 * Formato: cada línea contiene "x y"
 * Retorna el número de puntos leídos
 */
int leerDatosDesdeArchivo(const char *filename, double **x, double **y)
{
    FILE *file = fopen(filename, "r");
    if (file == NULL) {
        printf("[ERROR] No se pudo abrir el archivo '%s'\n", filename);
        return 0;
    }
    
    /* Contar líneas */
    int n = 0;
    char line[1024];
    while (fgets(line, sizeof(line), file)) {
        /* Ignorar líneas vacías o que comienzan con # */
        if (line[0] != '\n' && line[0] != '#') {
            n++;
        }
    }
    rewind(file);
    
    if (n == 0) {
        printf("[ERROR] El archivo está vacío o no tiene datos válidos\n");
        fclose(file);
        return 0;
    }
    
    /* Asignar memoria */
    *x = (double *)malloc(n * sizeof(double));
    *y = (double *)malloc(n * sizeof(double));
    
    if (!(*x) || !(*y)) {
        printf("[ERROR] Error de asignación de memoria\n");
        free(*x);
        free(*y);
        fclose(file);
        return 0;
    }
    
    /* Leer datos */
    int i = 0;
    while (fgets(line, sizeof(line), file) && i < n) {
        if (line[0] != '\n' && line[0] != '#') {
            if (sscanf(line, "%lf %lf", &(*x)[i], &(*y)[i]) == 2) {
                i++;
            }
        }
    }
    
    fclose(file);
    
    if (i != n) {
        printf("[ADVERTENCIA] Se esperaban %d puntos pero se leyeron %d\n", n, i);
        return i;
    }
    
    return n;
}