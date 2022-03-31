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
                 vida_media,genotypeFreq,gen=0,parents=0):
        
        #quiza el "name"no sea necesario para nada
        self.spName = name
        self.spSize = size
        self.spPloidy = ploidy

        self.vida_media = vida_media
        self.genotypeFreq = genotypeFreq
        
        self.sex = self.sex()
        self.ide = 'g'+str(gen)+".ID-"+str(nom)
        self.age = 0

        self.parents = parents
        self.genotype = dict()
        
        self.firstGenotype()
    
        

        
    def sex(self):
        if randint(0,1)==0:
            return "Male"
        else:
            return 'Female'
        
    def edad(self):
        #el usuario debera pasar algun parametro para indicar la dist por edades
        pass
    
    def firstGenotype(self,homogeneous=False):
        genotype = {}
        if self.parents == 0:

            if homogeneous:
                genotype = {'gene_1':'Aa','gene_2':'Bb'}
            else:
                genotype['gene_1']= ''.join(random.choices(Individual.gen1List,
                                            weights=self.genotypeFreq['A'],k=1))
                genotype['gene_2']= ''.join(random.choices(Individual.gen2List,
                                            weights=self.genotypeFreq['B'],k=1))
            print('hola')
            self.genotype = genotype
        else:
            self.mating()
    
    #metodo dunder
    def __str__(self):
        return ("este individuo es {}, su sexo es {} y su genotipo es {}"
              .format(self.ide,self.sex,self.genotype))
    
    def printParents(self,):
        '''
        Metodo que printa los padres del primer individuo
        '''
        print(f'''su padre es {self.parents[0].ide}, 
        con genotipo {self.parents[0].genotype}\n su madre es {self.parents[1].ide}
        con genotipo {self.parents[1].genotype}''')
        
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
        
            self.genotype['gene_1'] = aleloA_p1+aleloA_p2
            self.genotype['gene_2'] = aleloB_p1+aleloB_p2

            

