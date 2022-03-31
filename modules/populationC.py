'''
Clase Poblacion, permite crear una nueva poblacion, hacerla evolucionar
y obtener resumen de sus caracteristicas.
'''

from individualC import Individual

import random
import itertools


#clase poblacion,atributos generales que heredara de los individuos
class Population:
    
    genotypeFreq = dict()

    
    def __init__(self,size = 100,name="Homo sapiens",ploidy = 2, vida_media=55,
                 R=0.1,mu = (1e-4,1e-4),freq={'A':(0.4,0.6),'B':(0.6,0.4)}):
                 
        self.name = name
        self.size = size
        self.ploidy = ploidy
        self.vida_media = vida_media

        #frecuencia genotipica inicial
        self.freq = freq
        self.mu = mu
        self.genotypeFreq = self.genotFreq()
        self.gen = 0
        
    def __str__(self):
        return ''.join([self.name])
    
    #genera indivividuos
    def generateIndividuals(self):
        #ejecuta la funcion para crear las variables globales que contengan los genotipos posibles
        self.indiv = [Individual(i,self.name,
                                 self.size,
                                 self.ploidy,
                                 self.vida_media,
                                 self.genotypeFreq) 
                      for i in range(self.size)]
        print("se han generado un total de {} individuos de la poblacion {}"
              .format(self.size,self.name))
        
    #printa individuos        
    def printIndiv(self,show=5,children=True):
        show = abs(show)
        listaAtrib = ['ide','sex','genotype']
        print(*listaAtrib,sep="\t")
        if children==True and hasattr(self,'childrenInd'):
            print("print chidren")
            objectList = self.childrenInd
        else:
            objectList = self.indiv
            
        for x in objectList:
            print (*[getattr(x,y) for y in listaAtrib],sep="\t")
            #contador inverso, si se han enseñado show elementos para la ejecucion
            show += -1
            if show == 0:
                break
    
    def getInfo(poblacion):
        '''Esta funcion recogera un resumen de la poblacion
        y almacenara en un diccionario de clave 'generacion X' y valor una tupla
        con la informacion de la frecuencia genotipica
        '''
        pass

                
    #funcion que calcula las frecuencias genotipicas a partir de las alelicas
    def genotFreq(self):
        if self.ploidy == 2:
            for key,lista in self.freq.items():
                self.genotypeFreq[key] = (lista[0]**2,
                                         lista[0]*lista[1]*2,
                                         lista[1]**2)
                self.genotypeFreq[key] = [x*self.size for x in self.genotypeFreq[key]]
        elif self.ploidy == 1:
            self.genotypeFreq = self.freq
        return self.genotypeFreq
        
    def getMeanAge(self):
        '''obtienes la edad media recorriendo la lista de individuos'''
        try:
            meanAge = 0
            for obj in self.indiv:
                meanAge += obj.age
            meanAge = round(meanAge/len(self.indiv),2)
            print("la edad media de la poblacion es: ",meanAge)
        except:
            print("No has inicializado la poblacion")
            
                
    #muestra informacion genotipo
    def printGenotype(self):
        counter = {'AA':0,'Aa':0,'aa':0,'BB':0,'Bb':0,'bb':0}      
        for individuo in self.indiv:
            for key in counter:
                if individuo.genotype['gene_1'] == key or individuo.genotype['gene_2']==key:
                    counter[key] += 1
        print(counter)



      
    def evolvePop(self,gens = 20,every=5):

        for veces in range(0,gens):
            #aumentamos la generacion
            self.gen += 1
            #hacemos que poblacion apunte a la lista padre
            poblacion = self.indiv
            #vaciamos la lista individuos
            self.childrenInd = []
            #va introduciendo nuevos individuos hasta llegar al size de la poblacion
            for x in range(self.size): 
                self.childrenInd.append(self.chooseMate(x, poblacion))

            #sobreescribimos la generacion padre por la hija
            self.indiv = self.childrenInd

            #cada x generaciones, printamos
            if self.gen % every == 0:
                #este print sera otro metodo para obtener un resumen de 
                #la poblacion
                self.printIndiv(show=5)
                self.printGenotype()
                #self.printParentIndividuals(3)
        
        
    def chooseMate(self,x,poblacion):
        # elige dos individuos de forma aleatoria
        ind1,ind2 = random.choices(poblacion,k=2) 
        #guardamos los dos individuos en la variable parents
        parents = ind1,ind2
        # genera las ''frecuencias'' alelicas de ambos individuos, SIN USAR
        #self.findFreqAlleles(ind1,ind2)
        # genera las ''frecuencias'' genotipicas, NO ES NECESARIO AQUI
        #self.genotFreq()
        # nuevo nombre que se le pasara al Individual

        Ind_Name = x
        # genera un nuevo individuo y lo devuelve al metodo evolvePop
        return Individual(Ind_Name,
                         self.name,
                         self.size,
                         self.ploidy,
                         self.vida_media,
                         self.genotypeFreq,
                         self.gen,
                         parents)
    
    # SIN USAR: la eleccion de genotipo se hace ahora en la clae individual
    def findFreqAlleles(self,ind1,ind2):
        '''vacia el diccionario freq y cuenta las A para ambos padres,las B para ambos...
        esto lo almacena como freq (frecuencias alelicas)
        esto hay que cambiarlo porque freq no puede variar con cada individuo'''
        self.freq = dict()
        
        sizeA = len(ind1.genotype['gene_1'])+len(ind2.genotype['gene_1'])
        freq_A = (ind1.genotype['gene_1'].count('A')+ind2.genotype['gene_1'].count('A'))/sizeA
        freq_a = ind1.genotype['gene_1'].count('a')+ind1.genotype['gene_1'].count('a')/sizeA
        self.freq['A']=(freq_A,freq_a)
        
        sizeB = len(ind1.genotype['gene_2'])+len(ind2.genotype['gene_2'])
        freq_B = (ind1.genotype['gene_2'].count('B')+ind2.genotype['gene_2'].count('B'))/sizeB
        freq_b = (ind1.genotype['gene_2'].count('b')+ind2.genotype['gene_2'].count('b'))/sizeB
        self.freq['B']=(freq_B,freq_b)

    def listIndividuals(self):
        pass

    def printParentIndividuals(self,id=0):
        print(self.indiv[id])
        self.indiv[id].printParents()


   
if __name__ == '__main__':
    #se crea una nueva poblacion donde se especifican caracteristicas generales de esta
    #size es el numero de individuos
    #name el nombre de la especie
    #ploidy es la ploidia de la poblacion (haploide=1,diploide=2)
    #vida media es la vida media
    #freq son las frecuencias alelicas en cada locus, es una tupla de diccionarios
    #mu es la tasa de mutacion de los alelos (de A a a y viceversa..)
    
    shark = Population(size=10,
                        name="Megadolon",
                        ploidy=2,
                        vida_media=23,
                        freq={'A':(0.4,0.6),'B':(0.6,0.4)})

    #se generan individuos en esa poblacion
    shark.generateIndividuals()

    #parametro opcional show, permite elegir cuantos elementos se muestran (por defecto se muestran 10)
    shark.printIndiv(show=10)

    #muestra la cantidad de individuos con 'AA','aa'...
    shark.printGenotype()

    shark.evolvePop()

    shark.printIndiv(show=10)

    #printa el individuo que se quiere estudiar y sus padres
    shark.printParentIndividuals(id=2)

    #este metodo si lo llama el usuario te dara la info de la ultima generacion
    # shark.getInfo(5)
