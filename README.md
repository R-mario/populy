# Populy

Populy es un paquete de Python que permite llevar a cabo una simulación de tipo *forward time*. 
El paquete consta de dos módulos principales llamados *population* e *individual*, dentro se encuentran sus respectivas clases que permiten llevar a cabo la creación de 0 de una población, la evolucion de ésta y la obtención de unos informes, gráficos y resultados de la evolución.


## Instalación
Por el momento, se puede acceder mediante la descarga en [github](https://github.com/R-mario/populy)
en el futuro se podrá instalar con el gestor de paquetes [pip](https://pip.pypa.io/en/stable/).
```cmd
pip install populy
```
## Uso
```python
from populy.population import Population

# crea un objeto de la clase Poblacion
pop = Population(size=1000, 
                ploidy=2)

# genera individuos en la poblacio'n
pop.generateIndividuals()

# hace evolucionar a la Poblacion
pop.evolvePop(gens=200)
```
Para una explicación más detallada consultar el notebook [example](example.ipynb)

## Licencia
[MIT](https://choosealicense.com/licenses/mit/)