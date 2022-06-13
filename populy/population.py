
#local imports
from .individual import Individual
from .functions import fitness,outer_product

import random
import itertools
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 
import re
import os
from IPython.display import clear_output




#clase poblacion,atributos generales que heredara de los individuos
class Population:


    def __init__(self,size = 100,name="Population",ploidy = 2, vida_media=55,
                 R=0.5,mu = (0,0),freq={'A':(0.5,0.5),'B':(0.5,0.5)},D=0,
                 fit=0,sex_system='XY',rnd=False):
        """Creates a new empty population object.
        
        Parameters:
            size (int): Population size. Defaults to 100.
            name (int): Population name. Defaults to 'Population'
            ploidy (int): Number of homologous chromosomes. Defaults to 2.
            R (float) : Recombination frequency [0,0.5] where 0.5 indicates
            statistic independence. Defaults to 0.5. 
            mu (tuple(float,float)): Mutation rate. Defaults to (1e-4,1e-4).
            freq (dict): loci (key) and allelic frequencies for each allele (values).
            D (float): initial linkage desequilibrium [0,0.5]. Defaults to 0. 
            fit (int,dict): fitness function applied, can take allele fitness value
            or genotype fitness value. E.g. {'A':0.8} or {'AABB':0.8}. Every other 
            value is set to 1. Defaults to 0 (no fitness function applied).
            sex_system (str) : sex determination sistem. Can be XY,ZW or X0. 
            Defaults to 'XY'.
            rnd (bool) : If set to True every other parameter is changed to
            a new random value (R,D,mu,freq). Defaults to False.
        """     
        self.name = name
        self.size = size
        self.ploidy = ploidy
        self.vida_media = vida_media
        self.d = D
        self.R = R
        self.steps = 1
        
        self.rnd = rnd
        
        #frecuencia alelica inicial
        self.freq = self.initialAlleles(freq)
        self.gen = 0
        
        # stops evolve if needed
        self.stopEv = False
        
        self.fit = fit
        
        self.sex_system = sex_system.upper()
        
        self.mu = self.lenMutFreq(mu)
        self.__checkRandom()
        
    def initialAlleles(self,freq):
        """
        Checks if dictionary has the correct format

        Parameters:
            freq (dict): contains allelic frequencies values

        Raises:
            ValueError: if sum of all freq is greater than 1

        Returns:
            dict: diccionario de frecuencias
        """
        for k,v in freq.items():
            if self.rnd:
                q = random.random()
                freq[k]=(q,1-q)
            elif sum(freq[k])>1:
                raise ValueError(f'{freq[k]} is greater than 1') 
            elif isinstance(v,int):
                freq[k]=(v,1-v)
        return freq
    
    def lenMutFreq(self,mu):
        freqSize = len(self.freq)
        muList = list(mu)
        while(freqSize > len(muList)): 
            muList.append(0)
        
        return tuple(muList)
        
    
    def __checkRandom(self):
        '''
        Check if user set random (rnd) to True, then changes R, D, 
        mu and fit to random values.
           
        '''
        rnd = self.rnd
        if isinstance(rnd,bool):
            if rnd:
                self.R = random.random()/2
                self.d = random.random()/2
                self.mu = tuple([random.random()/2 for x in range(len(self.freq))])
                self.fit = random.randint(0,3)
                
    def __str__(self):
        return ''.join([self.name])
    
    
    def initIndividuals(self,pop=0):
        '''
        Creates new list of individuals
        
        Parameters:
            pop : If a population is passed then will use
            that as a generation 0 population. Defaults to 0.
        '''
        if not pop:     
            self.individuals = [Individual(i,
                                    self.name,
                                    self.size,
                                    self.ploidy,
                                    self.vida_media,
                                    self.freq,
                                    self.d,
                                    self.R,
                                    self.mu,
                                    self.sex_system,
                                    self.gen,
                                    parents=0) 
                        for i in range(self.size)]
            print(f"Se han generado {self.size} individuos en base a los atributos de la población")
        elif isinstance(pop,Population):
            # la poblacion sera una parte o toda, dependiendo del tamaño original
            # ERROR: si se quiere crear una población de una población
            self.individuals = Population.getCurrentIndividuals(pop)
            print(f"La población se ha iniciado con {len(self.individuals)} individuos de la población introducida ")
        elif isinstance(pop,list) and isinstance(pop[0],Individual):
            self.individuals = pop
            print(f"La población se ha iniciado con {len(self.individuals)} individuos introducidos ")
        else:
            print(f"{type(pop)} no es válido")
            
        
        # se crean nuevas variables de la poblacion
        # print(self.freq,self.alleleFreq())
        freq_alelicas = self.alleleFreq()
        # frecuencia alelica acumulada = se añadiran valores durante la ev
        self.f_ale_acc = {k: [v] for k,v in freq_alelicas.items()}
        dictc = self.gameticFreq()
        # frecuencia gametica acumulada
        self.f_gam_acc = {k: [v] for k,v in dictc.items()}
        # frecuencia de mutacion acumulada
        self.f_mut_acc = [self.findMutated()]
        # frecuencia de sexos acumulada
        self.f_sex_acc = [self.sexFreq()]
        
           
    def printIndividuals(self,show=5,children=True):
        '''
        Shows information about the individuals. 
        
        Parameters:
            show (int): Indicates how many individuals are shown. Defaults
            to 5. 
            children (bool): If True, shows current generation, if false,
            shows previous generation (parent). Defaults to True.
        '''
        show = abs(show)
        listaAtrib = ['ide','sex','sex_chromosome','chromosome']
        print(*listaAtrib,sep="\t")
        if children==True and hasattr(self,'childrenInd'):
            print("print chidren")
            objectList = self.childrenInd
        else:
            objectList = self.individuals
            
        for x in objectList:
            print (*[getattr(x,y) for y in listaAtrib],sep="\t")
            # contador inverso, si se han enseñado show elementos para la ejecucion
            show += -1
            if show == 0:
                break
    
    def plotAll(self):
        '''
        Graphical representation with matplotlib.
        Allelic and gametic frequencies, sex frequency and number of mutations.
        '''  

        # creamos el indice (eje x) del dataFrame
 
        labels = ['gen.'+str(x) for x in range(0,self.gen+1,self.steps)]
        
        # DataFrame de frecuencias alelicas acumuladas
        al_df = pd.DataFrame(self.f_ale_acc,index=labels)
        # DataFrame de frecuencias gameticas acumuladas
        gam_df = pd.DataFrame(self.f_gam_acc,index=labels)
        
        # DataFrame de sexos
        sex_df = pd.DataFrame(self.f_sex_acc,index=labels)

        # Hacemos el grafico
        fig,ax = plt.subplots(2,2,figsize=(13,8))
        
        plt.style.context("ggplot")
        plt.suptitle(f"Population with {self.size} individuals",fontsize=18)
        caption=f"""Initial conditions: allelic f.={self.freq}, recombination f.={self.R}
        mutation f.={self.mu}"""
        plt.figtext(0.5, 0.01, caption, wrap=True, horizontalalignment='center', fontsize=10)
        
        # fig[0].title('Variacion de las frecuencias gameticas')
        with plt.style.context("ggplot"):
          ax[0,0].plot(gam_df)
          ax[0,0].set_title('gametic frequencies')
          ax[0,0].legend(gam_df.columns)
          ax[0,0].set_ylim(0,1)
          ax[0,0].set_ylabel('p_gamete')
          ax[0,0].set_ylim([0,1])
          # fig[1].title('Variacion de las frecuencias alelicas')
          ax[0,1].set_title('Major allele frequencies')
          ax[0,1].plot(al_df)
          ax[0,1].legend(al_df.columns)
          ax[0,1].set_ylim(0,1)
          ax[0,1].set_ylabel('p_allele')
          ax[0,1].set_ylim([0,1])
          
          ax[1,0].plot(sex_df)
          ax[1,0].set_title('Sex frecuencies')
          ax[1,0].set_ylim(0.3,0.7)
          ax[1,0].legend(['Female','Male'])
          ax[1,0].set_ylabel('f(sex)')
          ax[1,0].set_ylim([0,1])
          
          ax[1,1].bar(height=self.f_mut_acc[1:],x=labels[1:])
          ax[1,1].set_title('Number of mutations')
          maxMutLim = max(self.f_mut_acc)
          ax[1,1].set_ylim(0,maxMutLim+1)
          ax[1,1].plot(self.f_mut_acc[1:], color='red', linewidth=2)
          ax[1,1].set_ylabel('mutations')
        
        new_steps = int(len(labels)/5) if len(labels)>8 else 1
        plt.setp(ax, xticks=range(0,len(labels),new_steps), xticklabels=labels[::new_steps])
        plt.show()

        
    def getDataFrame(self,which='mutantes'):
        '''
        Generates a pandas dataframe.
        
        Parameters:
            which (str) : which dataframe will be returned.
            Options = gametic,allelic,sex frequencies or number 
            of mutations
        
        Returns:
            (pd.DataFrame): Dataframe.
        '''
        labels = ['gen.'+str(x) for x in range(0,self.gen+1,self.steps)]
        
        if isinstance(which, str):
            if re.match('(.+?)?gamet(ica|o|e)s?',which):
                data = self.f_gam_acc
                Summary=pd.DataFrame(data,index=labels)
                Summary.columns = [f'p({i})' for i in Summary.columns]
            elif re.match('(.+?)?all?el(ica|o|e)s?',which):
                data = self.f_ale_acc
                Summary = pd.DataFrame(data,index=labels)
                Summary.columns = [f'p({i})' for i in Summary.columns]
            elif re.match('(.+?)?sexo?s?',which):
                data = self.f_sex_acc
                Summary =pd.DataFrame(data,index=labels,columns=['Female','Male'])
            elif re.match('(.+?)?mut(.+?)',which):
                data = self.f_mut_acc
                Summary = pd.DataFrame(data,index=labels)
            else:
                raise ValueError(f'Unknown {which}')
        elif isinstance(which,list):
            data = which
            Summary = pd.DataFrame(data,index=labels)
        else:
            raise TypeError(f'Unknown {which}') 

        return Summary
                                  
    
    def gameticFreq(self):
        '''
        Computes number of different gametes in the population.
        '''
        # diccionario tipo {'AB': 0,'Ab':0,...}
        obsGamf = outer_product(self.freq)
        obsGamf = {k:0 for k in obsGamf.keys()}
        # cuenta las ocurrencias en la poblacion de los distintos genotipos gameticos
        for individuo in self.individuals:
            for key in obsGamf:
                if individuo.chromosome['c1'] == key:
                    obsGamf[key] += 1
                if individuo.chromosome['c2']==key:
                    obsGamf[key] += 1
                    
        return {k: v / (2*len(self.individuals)) for k, v in obsGamf.items()}
    
    def setGameticFreq(self,gamf):
        self.gameticFrequencies = gamf
    
    
    def alleleFreq(self):
        '''
        Obtains frequency of major alleles from gametic frequency for current population.
        
        Returns:
            dict (str:int): Keys = locus(A,B,...), Values=Count
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
    
    def setAllelicFreq(self,alef):
        self.allelicFrequencies = alef
        
    def sexFreq(self):
        """Obtains sex frequency for current population.

        Returns:
            dict: Keys = Male-Female, Values = frequencies
        """
        sex=[0,0]
        for individuo in self.individuals:
            if individuo.getSex()=='Female':
                sex[0]+=1
            else:
                sex[1]+=1
        return [i/len(self.individuals) for i in sex]
    
            
    
    def freqGamAcumulada(self):
        '''
        Adds current gametic frequency to object variable f_gam_acc
        '''

        obsGamf = self.gameticFreq()
        self.setGameticFreq(obsGamf)
        # print(f'Generacion {self.gen}','frecuencia absoluta: ',obsGamf,sep='\n')
        # print(self.cum_gamF)

        # frecuencias gameticas acumuladas (durante las generaciones)
        for k in obsGamf:
            self.f_gam_acc[k].append(obsGamf[k])
    
    def freqAleAcumulada(self):
        '''
        Adds current allelic frequency to object variable f_ale_acc
        '''
        obsAleF = self.alleleFreq()
        self.setAllelicFreq(obsAleF)
        for k in obsAleF:
            self.f_ale_acc[k].append(obsAleF[k])
        
    def sexAcumulada(self):
        '''
        Adds current sex frequency to object variable f_sex_acc
        '''
        self.f_sex_acc.append(self.sexFreq())

    def mutAcumulada(self):
        '''
        Adds current mutation numbers to object variable f_mut_acc
        '''
        self.f_mut_acc.append(self.findMutated()) 
        
    def getObsGenotype(self,locus):
        """
        Obtains observed genotypes for current population.
        
        Parameters:
            locus (str): Locus of interest.
            
        Returns:
            dict: Keys = genotype, Values = number of individuals with that genotype.
        """
        try:
            
            genotype = {f'{locus}{locus}':0,f'{locus}{locus.lower()}':0,f'{locus.lower()}{locus.lower()}':0}
            
            for ind in self.individuals:
                for x in range(0,len(ind.genotype)):
                    if ind.genotype[x] in genotype.keys():
                        genotype[ind.genotype[x]] +=1
                        
            return genotype
        except:
            print(f'Wrong locus, it should be capitalize and one of {self.freq.keys()}')
        
    def getExpGenotype(self,locus,what_pop='initial'):
        '''Returns the expected genotype based on allele frequency and population size.
        The main use is for testing Hardy-Weinberg equilibrium.
        
        Parameters:
            locus (str): Locus of interest.
            what_pop (str): What population we are using to "expect" the frequencies. Values are "initial",
            which means that the gen.0 allele frecuencies will be to compute the value instead of current ones, 
            or "current". Defaults to "initial".
        
        Returns:
            dict: Keys = genotype, Values = number of expected individuals with that genotype.
        '''
        if what_pop == 'current':
            alF = self.alleleFreq()
        else:
            alF = {k:v[0] for k,v in self.freq.items()}
        
        try:
            p = alF[locus]
            q = 1 - p
            expAlFreq = {f'{locus}{locus}':p*p,f'{locus}{locus.lower()}':2*p*q,f'{locus.lower()}{locus.lower()}':q*q}
        except:
            print(f'Wrong locus, it should be capitalize and one of {list(self.freq.keys())}')
        
        #returns total expected count for given size
        return {k:round(v*self.size,2) for k,v in expAlFreq.items()}
        
                        
    def evolvePop(self,gens = 50,every=10,ignoreSex=False,printInfo=False,fit=None):
        """
        Evolves the population 

        Parameters:
            gens (int, optional): Number of generations. Defaults to 50.
            every (int, optional): How often it will get information. Defaults to 10.
            ignoreSex (bool, optional): Ignore sex when generating new children
            (only set to True when posize is really small). Defaults to False.
            printInfo (bool, optional): Show process info like some individuals. Defaults to False.
        """
        
        self.steps = every
        if fit is not None:
            self.fit = fit
        for veces in range(0,gens):
            # si hay que parar la evolucion por algun motivo, sale del bucle
            if self.stopEv:
                print(f'Se ha detenido la evolucion en la generacion {self.gen}')
                break
            #aumentamos la generacion
            self.gen += 1
            #hacemos que poblacion apunte a la lista padre
            currentPop = self.individuals
            #vaciamos la lista individuos
            self.childrenInd = []

            # introduce nuevos individuos hasta llegar al size de la poblacion
            x = 0
            while len(self.childrenInd) < self.size:
                child = self.__chooseMate(x, currentPop, ignoreSex)
                # aplicamos una funcion fitness
                if fitness(self.fit,child.genotype) == True:
                    self.childrenInd.append(child)
                    x+=1

            #sobreescribimos la generacion padre por la hija
            if self.stopEv == False:
                self.individuals = self.childrenInd

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
                
                completed = ((veces+1)/gens)*100
                if completed < 100:
                    #os.system('cls') error al ejecutar en rmarkdown
                    #clear_output(wait=True)
                    print(f"{round(completed,1)}% completado...")
        else:
            #clear_output(wait=True)
            #os.system('cls')
            # self.getInfo()
            print("¡Evolucion completada!")
                
        
        
    def __chooseMate(self,x,currentPop,ignoreSex):
        """
        Choose two parents and generate new children.

        Args:
            x (int): used to name the individual
            currentPop (list): list of individual instances.
            ignoreSex (bool): if sex will be ignored.

        Returns:
            Individual : new individual (children)
        """
        # elige dos individuos de forma aleatoria
        ind1,ind2 = random.choices(currentPop,k=2)
        count = 0
        # si son del mismo sexo vuelve a elegir, se establece un limite al bucle por si es infinito
        # Esto puede pasar cuando solo hayan machos o hembras en una poblacion pequeña
        while ind1.sex_chromosome == ind2.sex_chromosome and count < 5*self.size and ignoreSex==False:
            ind1,ind2 = random.choices(currentPop,k=2)
            # comprueba que sean de sexos distintos
            count +=1
        # si siguen siendo del mismo sexo, entonces hay que parar
        if ind1.sex_chromosome == ind2.sex_chromosome and ignoreSex==False:
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
                         self.sex_system,
                         self.gen,
                         parents)
    
    def getInfo(self):
        '''
        Call other methods which obtain information about the population
        '''
        self.freqGamAcumulada()
        self.freqAleAcumulada()
        self.sexAcumulada()
        self.mutAcumulada()

    def printSummary(self):
        """
        Get summary about current population.
        
        Will be obsolete.
        """
        tam = len(self.individuals)

        sex = {'Male':0,'Female':0}
        for x in range(tam):
            sexo = self.individuals[x].sex()
            if sexo == 'Male':
                sex['Male'] = sex['Male'] + 1
            else:
                sex['Female'] =sex['Female']+ 1

        print(f'Hay {len(self.individuals)} individuos\n{sex} son machos\t',
                f'{sex} son hembras \n\n el desequilibrio de ligamiento (LD) =',
                f'{self.d} \n frecuencia de recombinacion = {self.R} ',
                f' la generacion es {self.gen} las frecuencias gameticas', 
                f'hasta esta generacion son {self.f_gam_acc}')

    def printParentIndividuals(self,id=0):
        """Shows an individual an its parents

        Parameters:
            id (int, optional): which individual will be shown. Defaults to 0.
        """
        print(self.individuals[id])
        self.individuals[id].printParents()

    
    def findMutated(self,show=0):
        """Find mutated individuals

        Parameters:
            show (int, optional): If this search will be shown.
            Defaults to 0.

        Returns:
            int: number of individuals which were mutated.
        """
        mutated = 0
        advMutated = {k:0 for k in self.freq.values()}
        for individuo in self.individuals:
            if individuo.isMutated:    
                mutated += 1
                if show > mutated:
                    print("¡Un individuo ha mutado!, ha ocurrido",
                          f"en la generación {self.gen}",
                          " y se trata de:")
                    print(individuo)
            # if individuo.adMutated:
            #     for i in individuo.adMutated:
            #         advMutated[i] += 1
                    
        return mutated

    def info(pop):
        info = {'tamaño':pop.size,
                'ploidía':pop.ploidy, 
                'frecuencias alelicas iniciales':pop.freq, 
                'desequilibrio de ligamiento':pop.d, 
                'frecuencia de recombinacion':pop.R, 
                'tasa de mutaciones':pop.mu, 
                'generación actual':pop.gen,
                'sistema de determinación del sexo': pop.sex_system,
                'tipo de selección': pop.fit}
        if hasattr(pop,'indiv'):
            info['frecuencias alélicas actuales'] = pop.alleleFreq()
            info['número de individuos'] = len(pop.individuals)
        
        stringInfo = '\n'.join([f'{key}: {value}' for key, value in info.items()])
        print(stringInfo)
        
    def getCurrentIndividuals(self,howMany=None,shuffle=False):
        '''Returns a list of current individuals
        
        Parameters:
            howMany (int, optional): how many individuals will be returned.
            shuffle (bool, optional): if the individuals will be shuffled.
        
        Returns:
            list: list of current individuals
        '''
        individuos = self.individuals.copy()
        
        if shuffle:
            random.shuffle(individuos)
 
        if howMany is not None and howMany < self.size:
            individuos = individuos[:howMany]

        return individuos
    
    def fixedLoci(self):
        '''Returns which loci, if any, are fixed
        
        Returns:
            list: list of fixed loci
        '''
        loci = list()
        for k,v in self.allelicFrequencies.items():
            if v==0:
                loci.append(k.lower())
            elif v==1:
                loci.append(k)
            else:
                loci.append(0)
            
        
        return loci
    
    #TODO: implementar esta funcion
    def __iter__(self):
        pass
    def __next__(self):
        pass
                

   
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
    
    pop = Population(size=100,
                        name="Megadolon",
                        ploidy=2,
                        vida_media=23,
                        freq={'A':(0.4,0.6),'B':(0.6,0.4)},
                        D = 0.1,
                        R=0.5,
                        mu =(0.1,0.1),
                        fit = {'aabb':0.2})

    # se generan individuos en esa poblacion
    pop.initIndividuals()


    # parametro opcional show, permite elegir cuantos elementos se muestran (por defecto se muestran 10)
    pop.printIndiv(show=5)

    # muestra la cantidad de individuos con 'AA','aa'...
    # shark.printSummary()

    pop.evolvePop(gens=200,every=10,printInfo=False,ignoreSex=False)

    # printa el individuo que se quiere estudiar y sus padres
    # pop.printParentIndividuals(id=2)
    df = pop.getDataFrame('gametos')
    print(df.head())
    # obtiene un resumen del cambio en la frecuencia alelica
    pop.plotAll()
    



