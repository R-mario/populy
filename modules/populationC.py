'''
Clase Poblacion, permite crear una nueva poblacion, hacerla evolucionar
y obtener resumen de sus caracteristicas.
'''

from individualC import Individual

import random
import itertools
import pandas as pd
import matplotlib.pyplot as plt


#clase poblacion,atributos generales que heredara de los individuos
class Population:
    
    genotypeFreq = dict()

    
    def __init__(self,size = 100,name="Homo sapiens",ploidy = 2, vida_media=55,
                 R=0.1,mu = (1e-4,1e-4),freq={'A':(0.4,0.6),'B':(0.6,0.4)},D=0.1):
                 
        self.name = name
        self.size = size
        self.ploidy = ploidy
        self.vida_media = vida_media
        self.d = D
        self.R = R

        #frecuencia genotipica inicial
        self.freq = freq
        self.mu = mu
        self.genotypeFreq = self.genotFreq()
        self.gen = 0

        # igual a freq pero es una lista (se usara para printar)
        self.alleleFreqs = {k: [v[0]] for k,v in freq.items()}
        
        
        
    def __str__(self):
        return ''.join([self.name])
    
    # genera indivividuos
    def generateIndividuals(self):
        # ejecuta la funcion para crear las variables globales que contengan los genotipos posibles
        self.indiv = [Individual(i,
                                self.name,
                                 self.size,
                                 self.ploidy,
                                 self.vida_media,
                                 self.genotypeFreq,
                                 self.freq,
                                 self.d,
                                 self.R,
                                 self.gen) 
                      for i in range(self.size)]
        print("se han generado un total de {} individuos de la poblacion {}"
              .format(self.size,self.name))
        
        # se crean nuevas variables de la poblacion
        self.cum_gamF = self.gameticFreq()
        
    # printa individuos        
    def printIndiv(self,show=5,children=True):
        show = abs(show)
        listaAtrib = ['ide','sex','chromosome']
        print(*listaAtrib,sep="\t")
        if children==True and hasattr(self,'childrenInd'):
            print("print chidren")
            objectList = self.childrenInd
        else:
            objectList = self.indiv
            
        for x in objectList:
            print (*[getattr(x,y) for y in listaAtrib],sep="\t")
            # contador inverso, si se han enseñado show elementos para la ejecucion
            show += -1
            if show == 0:
                break
    
    def plotInfo(self):
        '''Esta funcion recogera un resumen del cambio en las frecuencias
        alelicas a lo largo de las generaciones, tambien obtiene un grafico del
        cambio. Se utilizan las librerias pandas y matplotlib
        '''
        df = pd.DataFrame(self.cum_gamF,index=[0])
        print(df)

        df.plot()
        plt.show()
        

                
    # funcion que calcula las frecuencias genotipicas a partir de las alelicas
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
            
                
    # calcula el numero de gametos distintos en la poblacion
    def gameticFreq(self):
        # frecuencia gametica observada
        obsGamf = {'AB':0,'Ab':0,'aB':0,'ab':0}    

        # cuenta las ocurrencias en la poblacion de los distintos genotipos  
        for individuo in self.indiv:
            for key in obsGamf:
                if individuo.chromosome['c1'] == key:
                    obsGamf[key] += 1
                if individuo.chromosome['c2']==key:
                    obsGamf[key] += 1
        
        return obsGamf
    
    def freqGamAcumulada(self):

        obsGamf = self.gameticFreq()

        print(f'Generacion {self.gen}','frecuencia absoluta: ',obsGamf,sep='\n')

        # frecuencias gameticas acumuladas (durante las generaciones)
        for k in obsGamf:
            self.cum_gamF[k] = self.cum_gamF[k].append(obsGamf[k])

        print(self.cum_gamF)
        # print('frecuencia alelica: ',self.alleleFreqs,sep='\n')



      
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
                # shark.printParentIndividuals(id=2)
                # self.printSummary()
                #self.printParentIndividuals(3)
        
        
    def chooseMate(self,x,poblacion):
        # elige dos individuos de forma aleatoria
        while True:
            ind1,ind2 = random.choices(poblacion,k=2)
            # comprueba que sean de sexos distintos
            if ind1.sex != ind2.sex:
                break 
        #guardamos los dos individuos en la variable parents
        parents = ind1,ind2
        # nuevo nombre que se le pasara al Individual

        Ind_Name = x
        # genera un nuevo individuo y lo devuelve al metodo evolvePop
        return Individual(Ind_Name,
                         self.name,
                         self.size,
                         self.ploidy,
                         self.vida_media,
                         self.genotypeFreq,
                         self.freq,
                         self.d,
                         self.R,
                         self.gen,
                         parents)

    def printSummary(self):
        tam = len(self.indiv)

        sex = {'Male':0,'Female':0}
        for x in range(tam):
            sexo = self.indiv[x].sex
            if sexo == 'Male':
                sex['Male'] = sex['Male'] + 1
            else:
                sex['Female'] =sex['Female']+ 1

        print(f'Hay {len(self.indiv)} individuos\n{sex} son machos\t',
                f'{sex} son hembras \n\n el desequilibrio de ligamiento (LD) =',
                f'{self.d} \n frecuencia de recombinacion = {self.R} ',
                f' la generacion es {self.gen} las frecuencias gameticas', 
                f'hasta esta generacion son {self.cum_gamF}')

    def printParentIndividuals(self,id=0):
        print(self.indiv[id])
        self.indiv[id].printParents()
    

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


   
if __name__ == '__main__':
    # se crea una nueva poblacion donde se especifican caracteristicas generales de esta
    # size es el numero de individuos
    # name el nombre de la especie
    # ploidy es la ploidia de la poblacion (haploide=1,diploide=2)
    # vida media es la vida media
    # freq son las frecuencias alelicas en cada locus, es una tupla de diccionarios
    # D es el desequilibrio de ligamiento de AB
    # R es la tasa de recombinacion
    # mu es la tasa de mutacion de los alelos (de A a a y viceversa..)
    
    shark = Population(size=5,
                        name="Megadolon",
                        ploidy=2,
                        vida_media=23,
                        freq={'A':(0.4,0.6),'B':(0.6,0.4)},
                        D = 0.1,
                        R=0)

    #se generan individuos en esa poblacion
    shark.generateIndividuals()


    #parametro opcional show, permite elegir cuantos elementos se muestran (por defecto se muestran 10)
    shark.printIndiv(show=5)

    #muestra la cantidad de individuos con 'AA','aa'...
    shark.printSummary()

    shark.evolvePop(gens=100,every=10)

    shark.printIndiv(show=5)

    #printa el individuo que se quiere estudiar y sus padres
    shark.printParentIndividuals(id=2)
    # obtiene un resumen del cambio en la frecuencia alelica
    shark.plotInfo()

