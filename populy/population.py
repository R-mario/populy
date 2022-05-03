'''
Clase Poblacion, permite crear una nueva poblacion, hacerla evolucionar
y obtener resumen de sus caracteristicas.
'''

from cProfile import label
from individual import Individual

import random
import itertools
import pandas as pd
import matplotlib.pyplot as plt
from functions import fitness,outer_product2


#clase poblacion,atributos generales que heredara de los individuos
class Population:
    
    genotypeFreq = dict()

    
    def __init__(self,size = 100,name="Population",ploidy = 2, vida_media=55,
                 R=0.1,mu = (1e-4,1e-4),freq={'A':(0.5,0.5),'B':(0.5,0.5)},D=0.1,
                 fit=0):
                 
        self.name = name
        self.size = size
        self.ploidy = ploidy
        self.vida_media = vida_media
        self.d = D
        self.R = R
        self.steps = 0

        #frecuencia genotipica inicial
        self.freq = self.allelicFreq(freq)
        self.mu = mu
        self.genotypeFreq = self.genotFreq()
        self.gen = 0
        self.mu = mu

        # igual a freq pero es una lista (se usara para printar)
        self.alleleFreqs = {k: [v[0]] for k,v in freq.items()}
        
        # variable booleana interrumpe la evolucion
        self.stopEv = False
        
        self.fit = fit
        
    def allelicFreq(self,freq):
        """Comprueba que el diccionario pasado tenga el formato correcto

        Args:
            freq (dict): contiene los valores de frecuencias alelicas

        Raises:
            ValueError: error si se pasan valores mayores que 1

        Returns:
            dict: diccionario de frecuencias
        """
        for k,v in freq.items():
            if type(v)==int:
                freq[k]=(v,1-v)
            if sum(freq[k])>1:
                raise ValueError('suma mayor que 1') 
        return freq
    
    
    def __str__(self):
        return ''.join([self.name])
    
    
    def generateIndividuals(self):
        '''
        Crea una lista de individuos
        '''       
        self.indiv = [Individual(i,
                                self.name,
                                 self.size,
                                 self.ploidy,
                                 self.vida_media,
                                 self.freq,
                                 self.d,
                                 self.R,
                                 self.mu,
                                 self.gen) 
                      for i in range(self.size)]
        print("se han generado un total de {} individuos de la poblacion"
              .format(self.size))
        
        # se crean nuevas variables de la poblacion
        dictc = self.gameticFreq()
        self.cum_gamF = {k: [v] for k,v in dictc.items()}
        print(self.cum_gamF)
        
           
    def printIndiv(self,show=5,children=True):
        '''
        Muestra por consola informacion de los primeros individuos
        de la poblacion 
        '''
        show = abs(show)
        listaAtrib = ['ide','sex','chromosome','isMutated']
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
    
    def plotInfo(self,what='Alleles',printInfo=False):
        '''
        Permite mostrar cambio de frecuencias geneticas a lo 
        largo de las generaciones.
        Se le puede indicar si se quiere mostrar 'alleles' o 'gametes'
        '''  

        if what=='Gametes':
            data = self.cum_gamF
        else:
            data = self.alleleFreqs
        # creamos el dataFrame
        index_name = ['gen.'+str(x) for x in range(0,self.gen+1,self.steps)]

        df = pd.DataFrame(data,index=index_name)

        al_df = pd.DataFrame(self.alleleFreqs,index=index_name)
        gam_df = pd.DataFrame(self.cum_gamF,index=index_name)
        if printInfo:
            print(df)
            print(al_df)

        # Hacemos el grafico
        fig,ax = plt.subplots(2,sharex=True, sharey=True)
        # fig[0].title('Variacion de las frecuencias gameticas')
        ax[0].plot(gam_df)
        ax[0].set_title('frecuencias gameticas')
        ax[0].legend(gam_df.columns)
        # fig[1].title('Variacion de las frecuencias alelicas')
        ax[1].set_title('frecuencias alelicas')
        ax[1].plot(al_df)
        ax[1].legend(al_df.columns)
        plt.show()
        
            
                
    
    def genotFreq(self):
        '''
        Calcula las frecuencias genotipicas a partir de las alelicas
        '''
        genotypeFreq = dict()
        if self.ploidy == 2:    
            for key,lista in self.freq.items():
                genotypeFreq[key] = (lista[0]**2,
                                         lista[0]*lista[1]*2,
                                         lista[1]**2)
                genotypeFreq[key] = [x*self.size for x in genotypeFreq[key]]
        elif self.ploidy == 1:
            genotypeFreq = self.freq
        return genotypeFreq

    # SIN USAR    
    def getMeanAge(self):
        '''
        obtienes la edad media recorriendo la lista de individuos
        '''
        try:
            meanAge = 0
            for obj in self.indiv:
                meanAge += obj.age
            meanAge = round(meanAge/len(self.indiv),2)
            print("la edad media de la poblacion es: ",meanAge)
        except:
            print("No has inicializado la poblacion")
            
                
    
    def gameticFreq(self):
        '''
        calcula el numero de gametos distintos en la poblacion
        '''
        # diccionario tipo {'AB': 0,'Ab':0,...}
        obsGamf = outer_product2(self.freq)
        obsGamf = {k:0 for k in obsGamf.keys()}
        # cuenta las ocurrencias en la poblacion de los distintos genotipos  
        for individuo in self.indiv:
            for key in obsGamf:
                if individuo.chromosome['c1'] == key:
                    obsGamf[key] += 1
                if individuo.chromosome['c2']==key:
                    obsGamf[key] += 1
                    
        return {k: v / (2*len(self.indiv)) for k, v in obsGamf.items()}
    
    def alleleFreq(self):
        '''
        Obtiene la frecuencia alelica a partir de la gametica para la poblacion actual
        '''
        # frecuencia alelica observada
        obsAleF = {k:0 for k in self.freq.keys()}
        # frecuencia gametica observada
        obsGamf = self.gameticFreq()
        for x in self.freq.keys():
            # suma los valores de frecuencia alelica que contengan la letra
            obsAleF[x] = sum(obsGamf[y] for y in obsGamf.keys() if x in y)
            
            # obsAleF['A'] = obsGamf['AB']+obsGamf['Ab']
            # obsAleF['B'] = obsGamf['AB']+obsGamf['aB']
            # ...

        return obsAleF
        

    
    def freqGamAcumulada(self):
        '''
        Modifica la variable cum_gamF para meter nuevos valores de frecuencia gametica
        '''

        obsGamf = self.gameticFreq()

        # print(f'Generacion {self.gen}','frecuencia absoluta: ',obsGamf,sep='\n')
        # print(self.cum_gamF)

        # frecuencias gameticas acumuladas (durante las generaciones)
        for k in obsGamf:
            self.cum_gamF[k].append(obsGamf[k])
    
    def freqAleAcumulada(self):
        '''
        Modifica la variable alleleFreqs para meter nuevos valores de frecuencia alelica
        '''
        obsAleF = self.alleleFreq()
        for k in obsAleF:
            self.alleleFreqs[k].append(obsAleF[k])
        


      
    def evolvePop(self,gens = 20,every=5,ignoreSex=True,printInfo=False):
        self.steps = every
        for veces in range(0,gens):
            # si hay que parar la evolucion por algun motivo, sale del bucle
            if self.stopEv:
                print(f'Se ha detenido la evolucion en la generacion {self.gen}')
                break
            #aumentamos la generacion
            self.gen += 1
            #hacemos que poblacion apunte a la lista padre
            poblacion = self.indiv
            #vaciamos la lista individuos
            self.childrenInd = []

            # introduce nuevos individuos hasta llegar al size de la poblacion
            x = 0
            while len(self.childrenInd)<= self.size:
                child = self.chooseMate(x, poblacion, ignoreSex)
                # aplicamos una funcion fitness
                if fitness(self.fit,child.genotype) == True:
                    self.childrenInd.append(child)
                    x+=1

            #sobreescribimos la generacion padre por la hija
            if self.stopEv == False:
                self.indiv = self.childrenInd

            # cada x generaciones, printamos
            if self.gen % every == 0:
                # enseña por pantalla informacion de individuos 
                # de la poblacion en curso si el usuario quiere
                if printInfo:    
                    self.printIndiv(show=5)
                
                # obtiene informacion de la poblacion en curso
                self.getInfo()
                
                # encuentra cuantos individuos han sufrido una mutacion
                self.findMutated(show = 2 if printInfo else 0)
                
                completed = (self.gen/gens)*100
                if completed < 100:
                    print(f"{round(completed,1)}% completado...")
        else:
            print("¡Evolucion completada!")
                
        
        
    def chooseMate(self,x,poblacion,ignoreSex):
        # elige dos individuos de forma aleatoria
        ind1,ind2 = random.choices(poblacion,k=2)
        count = 0
        # si son del mismo sexo vuelve a elegir, se establece un limite al bucle por si es infinito
        # Esto puede pasar cuando solo hayan machos o hembras en una poblacion pequeña
        while ind1.sex == ind2.sex and count < 5*self.size and ignoreSex==False:
            ind1,ind2 = random.choices(poblacion,k=2)
            # comprueba que sean de sexos distintos
            count +=1
        # si siguen siendo del mismo sexo, entonces hay que parar
        if ind1.sex == ind2.sex and ignoreSex==False:
            self.stopEv = True
           
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
                         self.freq,
                         self.d,
                         self.R,
                         self.mu,
                         self.gen,
                         parents)
    
    def getInfo(self):
        '''
        Llama a otros metodos que obtienen estadisticos de la poblacion
        '''
        self.freqGamAcumulada()
        self.freqAleAcumulada()

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

    
    def findMutated(self,show=10):
    # ver si algun individuo de la poblacion esta mutado 
        mutated = 0
        for individuo in self.indiv:
            if individuo.isMutated:    
                mutated += 1
                if show > mutated:
                    print("¡Un individuo ha mutado!, ha ocurrido",
                          f"en la generación {self.gen}",
                          " y se trata de:")
                    print(individuo)
        return mutated



   
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
    
    shark = Population(size=500,
                        name="Megadolon",
                        ploidy=2,
                        vida_media=23,
                        freq={'A':(0.4,0.6),'B':(0.6,0.4)},
                        D = 0.1,
                        R=0,
                        mu =(0,0),
                        fit = 3)

    # se generan individuos en esa poblacion
    shark.generateIndividuals()


    # parametro opcional show, permite elegir cuantos elementos se muestran (por defecto se muestran 10)
    shark.printIndiv(show=5)

    # muestra la cantidad de individuos con 'AA','aa'...
    # shark.printSummary()

    shark.evolvePop(gens=55,every=10,printInfo=False)

    shark.printIndiv(show=5)

    # printa el individuo que se quiere estudiar y sus padres
    shark.printParentIndividuals(id=2)
    # obtiene un resumen del cambio en la frecuencia alelica
    shark.plotInfo('Gametes')



