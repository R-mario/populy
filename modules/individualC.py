'''
Clase Individual, solo es llamada por la clase poblacion y contiene
informacion sobre un individuo concreto
'''
import itertools
from random import randint
import random

class Individual():

    def __init__(self,nom,name,size,ploidy,gen1List,gen2List,
                 vida_media,geneticPool):
        
        #quiza el "name"no sea necesario para nada
        self.spName = name
        self.spSize = size
        self.spPloidy = ploidy
        self.gene_1_List = gen1List
        self.gene_2_List = gen2List
        self.vida_media = vida_media
        self.geneticPool = geneticPool
        
        self.sex = self.sex()
        self.ide = "id-"+str(nom)
        self.age = 0
        self.genotype = self.firstGenotipo() #or self.generateGen()
        
    def sex(self):
        if randint(0,1)==0:
            return "Male"
        else:
            return 'Female'
        
    def edad(self):
        #el usuario debera pasar algun parametro para indicar la dist por edades
        pass
    
    def firstGenotipo(self,homogeneous=False):
        if homogeneous:
            return {'gen1':'Aa','gen2':'Bb'}
        else:
            return self.chooseGenotype()
        
    def chooseGenotype(self):
        dict = {}

        dict['gene_1']= ''.join(random.choices(self.gene_1_List,weights=self.geneticPool['A'],k=1))
        dict['gene_2']= ''.join(random.choices(self.gene_2_List,weights=self.geneticPool['B'],k=1))
        
        return dict
    
#en un principio esto no es necesario ya que no interesa conocer a un individuo concreto
    def __str__(self):
        return ("este individuo es {}, su sexo es {} y su genotipo es {}"
              .format(self.ide,self.sex,self.genotype))
    
    def childrenGen(self,list_a,list_b):
        pass
