from modules.populationC import Population
from modules.individualC import Individual

#se crea una nueva poblacion donde se especifican caracteristicas generales de esta
    #size es el numero de individuos
    #name el nombre de la especie
    #ploidy es la ploidia de la poblacion (haploide=1,diploide=2)
    #vida media es la vida media
    #freq son las frecuencias alelicas en cada locus, es una tupla de diccionarios
    #mu es la tasa de mutacion de los alelos (de A a a y viceversa..)
    
shark = Population(size=100,
                     name="Megadolon",
                     ploidy=2,
                     vida_media=23,
                     freq={'A':(0.4,0.6),'B':(0.6,0.4)})

#se generan individuos en esa poblacion
shark.generateIndividuals()

#parametro opcional show, permite elegir cuantos elementos se muestran (por defecto se muestran 10)
shark.printIndiv(show=10)

shark.evolvePop()
