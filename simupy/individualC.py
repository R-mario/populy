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
import itertools
from random import randint
import random

class Individual():

    gen1List = [''.join(x) for x in itertools.combinations_with_replacement('Aa',2)]
    gen2List = [''.join(x) for x in itertools.combinations_with_replacement('Bb',2)]


    def __init__(self,nom,name,size,ploidy,
                 vida_media,genotypeFreq,freq,d,R,mu,gen=0,
                 parents=0):
        
        # 
        self.spName = name
        self.spSize = size
        self.spPloidy = ploidy

        self.vida_media = vida_media
        self.genotypeFreq = genotypeFreq
        
        self.ide = 'g'+str(gen)+".ID-"+str(nom)
        self.age = 0
        self.parents = parents


        self.d = d
        self.alFreq = freq
        self.R = R
        self.mu = mu

        self.createIndividual()
        

        
    def createIndividual(self):
        '''
        Inicializa las variables que no se le pasan al inicializador
        '''
        self.sex = self.sex()
        self.genotype = dict()
        self.chromosome = dict()
        self.isMutated = False

        if self.parents:
            # self.mating()
            self.newMating()
            self.mutation()
        else:
            self.gameticFreq()


            
    def sex(self):
        if randint(0,1)==0:
            return "Male"
        else:
            return 'Female'
        
    def edad(self):
        #el usuario debera pasar algun parametro para indicar la dist por edades
        pass


    # Calcula la frecuencia gametica a partir de las frecuencias alelicas y D
    def gameticFreq(self):
        f = self.alFreq
        d = self.d
        fGametes = dict()
        fGametes['AB'] = f['A'][0]*f['B'][0]+d
        fGametes['Ab'] = f['A'][0]*f['B'][1]-d
        fGametes['aB'] = f['A'][1]*f['B'][0]-d
        fGametes['ab'] = f['A'][1]*f['B'][1]+d

        return self.chooseGametes(fGametes)

    # Elige un gameto segun su probabilidad (frecuencia)
    def chooseGametes(self,fGametes):
        chromosome = dict()
        gameto =list(fGametes.keys())
        pesos = list(fGametes.values())
        for i in range(1,self.spPloidy+1):
            chromosome['c'+str(i)]= ''.join(random.choices(gameto,
                                            weights=pesos,k=1))

        self.chromosome = chromosome


    
    # metodo dunder
    def __str__(self):
        return ("este individuo es {}, su sexo es {} y su genotipo es {}"
              .format(self.ide,self.sex,self.chromosome))
    
    def printParents(self,):
        '''
        Metodo que printa los padres del primer individuo
        '''
        print(f'''su padre es {self.parents[0].ide}, 
        con genotipo {self.parents[0].chromosome}\n su madre es {self.parents[1].ide}
        con genotipo {self.parents[1].chromosome}''')

    # calcula
    def newMating(self):
        r = self.R
        recomb = ((1-r)/2,(1-r)/2,r/2,r/2)
        for x in range(len(self.parents)):
            c1P = self.parents[x].chromosome['c1']
            c2P = self.parents[x].chromosome['c2']

            ch_P = [c1P,c2P,c1P[0]+c2P[1],c2P[0]+c1P[1]]
            # print('el padre tiene estas combinaciones: ',ch_P,recomb)
            self.chromosome['c'+str(x+1)]= random.choices(ch_P,weights=recomb,k=1)[0]
            # print('da lugar a este chromosoma: ',self.chromosome['c'+str(x+1)])
        # print(self.ide)
        # print('el resultado es: ',self.chromosome)
        
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

    # SIN USAR
    def mating(self):
        '''Obtiene un genotipo a partir de los padres'''
        p1_genotype = self.parents[0].genotype
        p2_genotype = self.parents[1].genotype

        if self.parents != 0:
            # para el alelo A
            if randint(0,1)==0:
                aleloA_p1 = p1_genotype['gene_1'][0]
                if randint(0,1)==0:
                    aleloA_p2 = p2_genotype['gene_1'][0]
                else:
                    aleloA_p2 = p2_genotype['gene_1'][1]
            else:
                aleloA_p1 = p1_genotype['gene_1'][1]
                if randint(0,1)==0:
                    aleloA_p2 = p2_genotype['gene_1'][0]
                else:
                    aleloA_p2 = p2_genotype['gene_1'][1]
            # para B
            if randint(0,1)==0:
                aleloB_p1 = p1_genotype['gene_2'][0]
                if randint(0,1)==0:
                    aleloB_p2 = p2_genotype['gene_2'][0]
                else:
                    aleloB_p2 = p2_genotype['gene_2'][1]
            else:
                aleloB_p1 = p1_genotype['gene_2'][1]
                if randint(0,1)==0:
                    aleloB_p2 = p2_genotype['gene_2'][0]
                else:
                    aleloB_p2 = p2_genotype['gene_2'][1]
        
            a = aleloA_p1+aleloA_p2
            b = aleloB_p1+aleloB_p2
            #se aplica esto para que todos los 'Aa' aparezcan asi en lugar de 'aA'
            self.genotype['gene_1'] = ''.join(sorted(a))
            self.genotype['gene_2'] = ''.join(sorted(b))

            

