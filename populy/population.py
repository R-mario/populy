'''
Clase Poblacion, permite crear una nueva poblacion, hacerla evolucionar
y obtener resumen de sus caracteristicas.
'''

from cProfile import label
from individual import Individual
from functions import fitness,outer_product

import random
import itertools
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 
from statistics import mean



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
        
        # frecuencia alelica acumulada = se añadiran valores durante la ev
        self.f_ale_acc = {k: [v[0]] for k,v in self.freq.items()}
        dictc = self.gameticFreq()
        # frecuencia gametica acumulada
        self.f_gam_acc = {k: [v] for k,v in dictc.items()}
        # frecuencia de mutacion acumulada
        self.f_mut_acc = [self.findMutated()]
        # frecuencia de sexos acumulada
        self.f_sex_acc = [self.sexFreq()]
        
           
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
    
    def plotAll(self,printInfo=False):
        '''
        Representa gráficamente todas las caracteristicas 
        más relevantes de la población 
        frecuencias alélicas, frecuencias gaméticas, sexo, nº mutantes y
        nº de recombinantes
        Args:
            printInfo (bool, optional): Mostrar o no mas informacion 
            a parte de las graficas. Defaults to False.
        '''  

        # creamos el indice (eje x) del dataFrame
 
        index_name = ['gen.'+str(x) for x in range(0,self.gen+1,self.steps)]
        
        # DataFrame de frecuencias alelicas acumuladas
        al_df = pd.DataFrame(self.f_ale_acc,index=index_name)
        # DataFrame de frecuencias gameticas acumuladas
        gam_df = pd.DataFrame(self.f_gam_acc,index=index_name)
        
        # DataFrame de sexos
        sex_df = pd.DataFrame(self.f_sex_acc,index=index_name)
        # # DataFrame de mutaciones 
        # mut_df = pd.DataFrame(self.f_mut_acc,index=index_name)
        if printInfo:
            print(gam_df)
            print(al_df)

        # Hacemos el grafico
        fig,ax = plt.subplots(2,2,figsize=(13,8))

        # fig[0].title('Variacion de las frecuencias gameticas')
        ax[0,0].plot(gam_df)
        ax[0,0].set_title('frecuencias gaméticas')
        ax[0,0].legend(gam_df.columns)
        ax[0,0].set_ylim(0,1)
        # fig[1].title('Variacion de las frecuencias alelicas')
        ax[0,1].set_title('frecuencias alélicas')
        ax[0,1].plot(al_df)
        ax[0,1].legend(al_df.columns)
        ax[0,1].set_ylim(0,1)
        
        ax[1,0].plot(sex_df)
        ax[1,0].set_title('frecuencia del sexo')
        ax[1,0].set_ylim(0.3,0.7)
        ax[1,0].legend(['Female','Male'])
        
        ax[1,1].bar(height=self.f_mut_acc[1:],x=index_name[1:])
        ax[1,1].set_title('Mutaciones')
        maxMutLim = max(self.f_mut_acc)
        ax[1,1].set_ylim(0,maxMutLim+1)
        # media
        ax[1,1].axhline(mean(self.f_mut_acc), color='green', linewidth=2)
        if len(index_name)>10:
            for a in np.ravel(ax): 
                a.set_xticks(index_name[1:-1:3])
        plt.show()
        
    def plotAlleles(freqs):
        '''
        La idea es que este metodo te permita representar el cambio 
        de las frecuencias alelicas si se llama por separado 
        pero si lo llama un metodo superior que se pueda meter en un
        subplot.
        Para eso debe ser un metodo estatico?
        '''
                         
    
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
        obsGamf = outer_product(self.freq)
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
        
    def sexFreq(self):
        sex=[0,0]
        for individuo in self.indiv:
            if individuo.sex=='Female':
                sex[0]+=1
            else:
                sex[1]+=1
        return [i/len(self.indiv) for i in sex]
    
            
    
    def freqGamAcumulada(self):
        '''
        Modifica la variable cum_gamF para meter nuevos valores de frecuencia gametica
        '''

        obsGamf = self.gameticFreq()

        # print(f'Generacion {self.gen}','frecuencia absoluta: ',obsGamf,sep='\n')
        # print(self.cum_gamF)

        # frecuencias gameticas acumuladas (durante las generaciones)
        for k in obsGamf:
            self.f_gam_acc[k].append(obsGamf[k])
    
    def freqAleAcumulada(self):
        '''
        Modifica la variable alleleFreqs para meter nuevos valores de frecuencia alelica
        '''
        obsAleF = self.alleleFreq()
        for k in obsAleF:
            self.f_ale_acc[k].append(obsAleF[k])
        
    def sexAcumulada(self):
        self.f_sex_acc.append(self.sexFreq())

    def mutAcumulada(self):
        self.f_mut_acc.append(self.findMutated()) 
        
         
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
        self.sexAcumulada()
        self.mutAcumulada()

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
                f'hasta esta generacion son {self.f_gam_acc}')

    def printParentIndividuals(self,id=0):
        print(self.indiv[id])
        self.indiv[id].printParents()

    
    def findMutated(self,show=0):
        """encuentra a los individuos que han mutado en la poblacion

        Args:
            show (int, optional): Cuántos individuos queremos que muestre.
            Defaults to 0.

        Returns:
            int: número de individuos que han sufrido una mutación
        """
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
    
    shark = Population(size=100,
                        name="Megadolon",
                        ploidy=2,
                        vida_media=23,
                        freq={'A':(0.4,0.6),'B':(0.6,0.4)},
                        D = 0.1,
                        R=0,
                        mu =(0.1,0.1),
                        fit = 3)

    # se generan individuos en esa poblacion
    shark.generateIndividuals()


    # parametro opcional show, permite elegir cuantos elementos se muestran (por defecto se muestran 10)
    shark.printIndiv(show=5)

    # muestra la cantidad de individuos con 'AA','aa'...
    # shark.printSummary()

    shark.evolvePop(gens=200,every=10,printInfo=False)

    shark.printIndiv(show=5)

    # printa el individuo que se quiere estudiar y sus padres
    shark.printParentIndividuals(id=2)
    # obtiene un resumen del cambio en la frecuencia alelica
    shark.plotAll()



