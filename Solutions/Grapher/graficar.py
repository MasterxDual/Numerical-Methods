#Paquetes
import numpy as np
import matplotlib.pyplot as plt

#Leo los datos
data = np.loadtxt('data.txt.dat') #leo el archivo
x = data[:,0] #primer columna
y = data[:,1] #segunda columna

#Grafico
plt.plot(x,y, marker='o', color='r') # ploteo el grafico
plt.xlim(min(x)-0.5, max(x)+0.5) #limites del eje x basados en los datos
plt.ylim(min(y)-0.5, max(y)+0.5) #limites del eje y basados en los datos
plt.title('plot') #titulo del grafico
plt.xlabel('x') #nombre del eje x
plt.ylabel('y') #nombre del eje y
plt.show()

"""
https://matplotlib.org/stable/api/markers_api.html ---> distintos tipos de puntos
https://matplotlib.org/stable/gallery/color/named_colors.html ----> colores
https://matplotlib.org/stable/api/_as_gen/matplotlib.lines.Line2D.html#matplotlib.lines.Line2D.set_linestyle --> tipo de linea
"""