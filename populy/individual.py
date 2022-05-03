'''
Clase Individual, solo es llamada por la clase poblacion y contiene
informacion sobre un individuo concreto
esto es:
str sexo: male or female
str ide: identificador que consta de generacion (g.X) + numero (ID-YYY)
int age: de momento sin usar
dict genotype: {'str': 'str'} siendo la clave A o B y los valores
                su genotipo para ese locus
'''
from calendar import c
import itertools
import numpy as np
from random import randint
import random
from functions import outer_product2

class Individual():


    def __init__(self,nom,name,size,ploidy,
                 vida_media,freq,d,R,mu,gen=0,
                 parents=0):
        
         
        self.spName = name
        self.spSize = size
        self.spPloidy = ploidy

        self.vida_media = vida_media
        
        self.ide = 'g'+str(gen)+".ID-"+str(nom)
        self.age = 0
        self.parents = parents


        self.d = d
        # diccionario con frecuencias alelicas
        self.alFreq = freq
        # frecuencia de recombinacion
        self.R = R
        # tasa de mutacion
        self.mu = mu

        self.createIndividual()
        

        
    def createIndividual(self):
        '''
        Inicializa las variables que no se le pasan al inicializador
        '''
        self.sex = self.sex()
        
        self.chromosome = dict()
        self.isMutated = False

        if self.parents:
            # self.mating()
            self.mating()
            self.mutation()
        else:
            self.chromosome = self.chooseGametes()
            self.gameticFreq()
        
        self.genotype = self.getGenotype()

            
    def sex(self):
        if randint(0,1)==0:
            return "Male"
        else:
            return 'Female'
        
    def getGenotype(self):
        
        genotype = dict()
        # print(self.chromosome)
        for i,letra in enumerate(self.chromosome['c1']):
            if self.spPloidy == 2:
                genotype[letra.upper()] = ''.join(sorted(self.chromosome['c1'][i] + self.chromosome['c2'][i]))
            else:
                genotype[letra.upper()] = self.chromosome['c1'][i]
        
        return genotype

    # Calcula la frecuencia gametica a partir de las frecuencias alelicas y D
    # esto quiza deberia estar en population !!!revisar
    def gameticFreq(self):

        f = self.alFreq
        d = self.d
        fGametes = dict()
        
        #producto externo (outer product)
        freqValues = list(f.values())
        inp = freqValues[0]
        finalDict = dict()
        if len(freqValues)>1:
            fGametes = outer_product2(f)
            # for x in range(1,len(freqValues)):
            #     final = outer_product(inp,freqValues[x],x,finalDict)
            #     inp = list(final.values())
            # fGametes = final
        else:
            word = str(list(f.keys())[0])
            keys = [word,word.lower()]
            
            values = list(freqValues[0])  
                 
            fGametes = dict(zip(keys,values))

        return fGametes

    # Elige un gameto segun su probabilidad (frecuencia)
    def chooseGametes(self):
        """Choose a gamete from gameticFreq keys
        by its given probability in gameticFreq values

        Returns:
           dict(str:str): two chromosomes and its genotype
        """
        chrom = dict()
        fGametes = self.gameticFreq()
        
        gameto =list(fGametes.keys())
        pesos = list(fGametes.values())
        
        for i in range(1,self.spPloidy+1):
            
            chrom['c'+str(i)]= ''.join(random.choices(gameto,
                                            weights=pesos,k=1))

        # self.chromosome = chrom
        return chrom

    
    # metodo dunder
    def __str__(self):
        return ("Este individuo es {}, su sexo es {} su genotipo es {}"
              .format(self.ide,self.sex,self.chromosome))
    
    def printParents(self,):
        '''
        Print individual parents
        '''
        
        print(f'''su padre es {self.parents[0].ide}, 
        con genotipo {self.parents[0].chromosome}\n su madre es {self.parents[1].ide}
        con genotipo {self.parents[1].chromosome}''')

    # calculates
    def mating(self):
        # print(self.chromosome,len(self.chromosome))
        if self.spPloidy > 1:
            r = self.R
            # tabla de recombinacion
            recomb = ((1-r)/2,(1-r)/2,r/2,r/2)
            for x in range(len(self.parents)):
                c1P = self.parents[x].chromosome['c1']
                c2P = self.parents[x].chromosome['c2']
                # para 2 locus (A,B)
                if len(self.alFreq) == 2:                 
                    ch_P = [c1P,c2P,c1P[0]+c2P[1],c2P[0]+c1P[1]]
                    self.chromosome['c'+str(x+1)]= random.choices(ch_P,weights=recomb,k=1)[0]
                # para 1 locus (A)
                elif len(self.alFreq)==1:
                    ch_P = [c1P,c2P]
                    self.chromosome['c'+str(x+1)]= random.choices(ch_P,weights=(0.5,0.5),k=1)[0]     
        else:  
            self.chromosome = {'c'+str(k+1):v for k in range(2) for v in random.choice(self.parents[k].chromosome.values())}

        
    def mutation(self):
        '''
        Provoca el cambio del alelo mayor al menor con una frecuencia mut
        '''
        muType = 'unidirectional'
        for k,v in self.chromosome.items():
            for i in range(len(v)):
                #comprueba si es el alelo mayor
                if v[i].isupper():
                    #si muta, cambia el alelo del cromosoma por el alelo menor
                    if random.random() < self.mu[i]:
                        self.chromosome[k] = self.chromosome[k].replace(v[i],v[i].lower())
                        self.isMutated = True

            

