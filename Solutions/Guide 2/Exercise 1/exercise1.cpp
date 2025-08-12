/*Sea g (x) = x2 + x − 4. Podemos utilizar iteración de punto fijo para hallar las soluciones de la
ecuación x = g (x)?. Por que?. Suponiendo una ecuación arbitraria que tiene un punto fijo P , Cual
es la ventaja de tener g'(P ) ≈ 0 en un proceso de iteración de punto fijo?.

Depende como sea hayamos tomado a la funcion dada, en este caso g(x) lo tomamos que es f(x) y luego conseguimos g(x) para hacer luego su derivada -> {
Segundo intento:
Despejamos x:
Tenemos que f(x) = x2 + x -4 = 0, luego x2 - 4 = -x, o equivalentemente x = 4 - x2, viendo que g(x) = 4 - x2. Entonces derivamos a g(x), teniendo
g'(x) = -2*x, luego fabs(g'(x)) = 2 * fabs(x) teniendo que fabs(g'(x)) > 1, llegando a la conclusion de que tampoco converge en ningun momento.

Tercer intento: 
Sumamos x a ambos miembros:
Nos queda x2 + 2*x - 4 = x, luego g(x) = x2 + 2*x -4, luego g'(x) = 2*x + 2, es decir g'(x) = 2*(x + 1), aplicando el valor absoluto nos queda
fabs(g'(x)) = 2 * fabs(x + 1), y esto quiere decir que fabs(g'(x)) > 1. Concluyendo que tampoco converge en ningun momento.
}

Si suponemos una ecuacion arbritaria que tiene un punto fijo P, la ventaja de tener g'(P) ≈ 0 quiere decir que g'(P) < 1, lo que nos quiere decir que la derivada 
evaluada en el punto fijo P nos dará un valor menor a 1, lo que nos dice que en cada iteracion nos hemos acercado cada vez mas al punto fijo P y ademas
hemos podido encontrar el punto fijo y que el programa pudo realizar
todos los calculos necesarios para llegar a la raiz de f(x) o al punto fijo de g(x), habiendo convergido la función. 
*/