#include <stdio.h>
#include <stdlib.h>

#define MATRIZ_TXT "matrix.txt"
#define N 100

int main(int argc, char const *argv[])
{
    FILE *file = fopen(MATRIZ_TXT, "w");
    if (file == NULL) {
        printf("The file could not be opened for writing.\n");
        return 1;
    }

    // Correct memory allocation for a 2D matrix
    double **A = (double **)malloc(N * sizeof(double *));
    double *b = (double *)malloc(N * sizeof(double));
    if (!A || !b)
    {
        printf("Memory allocation error\n");
        fclose(file);
        return 1;
    }

    for (size_t i = 0; i < N; i++)
    {
        A[i] = (double *)malloc(N * sizeof(double));
        if (!A[i]) {
            printf("Memory error in row %zu\n", i);
            // Here you would need to free the memory already allocated before exiting.
            fclose(file);
            return 1;
        }
        for (size_t j = 0; j < N; j++)
        {
            A[i][j] = 0;
        }
    }

    for (size_t i = 0; i < N; i++)
    {
        // Main Diagonal
        A[i][i] = -2;

        // Independent vector
        b[i] = 4;

        // Subdiagonal
        if (i > 0)
            A[i][i-1] = 1;

        // Superdiagonal
        if (i < N-1)
            A[i][i+1] = 1;
    }

    A[0][0] = 1;
    A[N-1][N-1] = 1;
    A[0][1] = 0;
    A[N-1][N-2] = 0;
    b[0] = 1;
    b[N-1] = 1;

    // Save the matrix and vector to the file
    printf("Saving %dx%d matrix to %s...\n", N, N, MATRIZ_TXT);
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            fprintf(file, "%.1lf ", A[i][j]);
        }
        // When you want to read the program, remember to remove the | and the tab.
        fprintf(file, "%.1lf\n", b[i]);
    }
    printf("Matrix saved successfully.\n");

    // Free memory
    fclose(file);
    for (size_t i = 0; i < N; i++) {
        free(A[i]);
    }
    free(A);
    free(b);

    return 0;
}