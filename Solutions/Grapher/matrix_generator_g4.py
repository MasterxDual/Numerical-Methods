# Generate data.dat for the band system in Problem 4
n = 50  # System size

"""
--> El índice i del for representa la fila actual que estás construyendo.
--> El [] dentro del row indica las columnas de la fila
--> Esto sucede porque primero definis una fila, luego se escribe esa fila, y vuelve a empezar el bucle for escribiendo otra fila debajo
--> La lista row representa el contenido de esa fila.
--> El write con \n asegura que cada fila esté en una línea diferente dentro del archivo. """
with open("data.dat", "w") as f:
    for i in range(n):
        row = [0]*n
        # Band coefficients
        if i-2 >= 0: 
            row[i-2] = 1
        if i-1 >= 0: 
            row[i-1] = -2
        row[i] = 12
        if i+1 < n: 
            row[i+1] = -2
        if i+2 < n: 
            row[i+2] = 1
        # Write row with its independent term (5)
        f.write(" ".join(f"{val:3d}" for val in row) + " 5\n")
