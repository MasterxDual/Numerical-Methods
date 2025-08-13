# Esto es un comentario
# La siguiente línea indica el compilador que se va a usar
CC = g++

# La siguiente línea indica las opciones del compilador
CFLAGS = -Wall

# La regla por defecto
all: method1

method1: method1.o
	$(CC) $(CFLAGS) method1.o -o method1

method1.o: 
	$(CC) $(CFLAGS) -c method1.cpp

run:
	./method1

# Regla para limpiar los archivos generados
clean:
	rm -f method1.o method1