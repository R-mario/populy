import unittest
import sys,os
from pathlib import Path

from populy import functions
# from populy.individual import Individual

try:
    import populy.individual as ind
    import populy.population as pop
except ModuleNotFoundError as e:
    path = Path()
    popuPath = os.path.join(os.path.dirname(__file__), 'populy')
    # print(popuPath)
    sys.path.append(popuPath)
    # print(sys.path)
    import populy.individual as ind
    import populy.population as pop


class Test_individual(unittest.TestCase):

    def test_individual_activation(self):
        # instanciamos
        ind1 = ind.Individual('nombre_ind','nombre_pob',1,2,0,freq={'A':(0.4,0.6),'B':(0.6,0.4)},
                          d=0,R=0.5,mu=(0.1,0.1))
        # comprobamos que devuelve un string
        self.assertTrue(ind1)
        
    def test_gametic(self):
        # instanciamos
        ind2 = ind.Individual('nombre_ind','nombre_pob',1,2,0,freq={'A':(0.4,0.6),'B':(0.6,0.4)},
                          d=0,R=0.5,mu=(0.1,0.1))
        # guardamos la frecuencia gametica
        gam_freq = ind2.gameticFreq()
        # comprobamos que es lo esperado
        self.assertDictEqual(gam_freq,{'AB': 0.24, 'Ab': 0.16000000000000003, 'aB': 0.36, 'ab': 0.24})
        
    def test_haploid(self):
        ind3 = ind.Individual('nombre_ind','nombre_pob',1,1,0,freq={'A':(0.4,0.6),'B':(0.6,0.4)},
                          d=0,R=0.5,mu=(0.1,0.1))
        #chromosomes 
        chr = ind3.chromosome
        self.assertEqual(list(chr.keys()),['c1'])
        genotype = ind3.genotype
        self.assertEqual(len(genotype['A']),1)
    
    def test_oneLocus(self):
        freq={'A':(0.4,0.6)}
        ind3 = ind.Individual('nombre_ind','nombre_pob',1,2,0,freq,
                          d=0,R=0.5,mu=(0.1,0.1))
        gamFreq=ind3.gameticFreq()
        self.assertEqual(gamFreq,{'A':0.4,'a':0.6})
    
    def test_mating_oneLocus(self):
        freq = {'A':(0.4,0.6)}
        # padres del individuo
        p1 = ind.Individual('nombre_ind','nombre_pob',1,2,0,freq,
                          d=0,R=0.5,mu=(0.1,0.1))
        p2 = ind.Individual('nombre_ind','nombre_pob',1,2,0,freq,
                          d=0,R=0.5,mu=(0.1,0.1))
        parents = [p1,p2]
        # individuo
        ind5 = ind.Individual('nombre_ind','nombre_pob',1,2,0,freq,
                          d=0,R=0.5,mu=(0.1,0.1),parents=parents)
        
        

class Test_population(unittest.TestCase):
    def test_pop_empty(self):
        """
        creates a population without giving parameters
        """
        pop1 = pop.Population()
        self.assertTrue(pop1)

    
    def test_pop_oneLocus(self):
        pop2 = pop.Population(freq={'A': (0.4,0.6),'B':(0.1,0.9)})
        pop2.generateIndividuals()
        pop2.evolvePop(printInfo=False)
        print(pop2)


class Test_functions(unittest.TestCase):
    def test_outer_three(self):
        f = {'A':0,'B':0,'C':0.3}
        newf = {k:(v,1-v) for k,v in f.items()}
        finalD = functions.outer_product2(newf)
        
        self.assertEqual(len(finalD),2**len(f))
    
    def test_outer_one(self):
        f = {'A':1,'B':1,'C':0.3}
        newf = {k:(v,1-v) for k,v in f.items()}
        fValues = list(newf.values())
        i = fValues[0]
        finalD = dict()
        finalD = functions.outer_product(i,fValues[1],1,finalD)
        # print('oldcic',finalD)

if __name__ == '__main__':
    unittest.main()
    