import unittest
import sys,os
from pathlib import Path

from populy import functions
# from populy.individual import Individual

try:
    import populy.individual as pop
except ModuleNotFoundError as e:
    path = Path()
    popuPath = os.path.join(os.path.dirname(__file__), 'populy')
    # print(popuPath)
    sys.path.append(popuPath)
    # print(sys.path)
    import populy.individual as pop


class individualTest(unittest.TestCase):

    def test_individual_activation(self):
        # instanciamos
        ind1 = pop.Individual('nombre_ind','nombre_pob',1,2,0,freq={'A':(0.4,0.6),'B':(0.6,0.4)},
                          d=0,R=0.5,mu=(0.1,0.1))
        # comprobamos que devuelve un string
        self.assertTrue(ind1)
    def test_gametic(self):
        ind2 = pop.Individual('nombre_ind','nombre_pob',1,2,0,freq={'A':(0.4,0.6),'B':(0.6,0.4)},
                          d=0,R=0.5,mu=(0.1,0.1))
        gam_freq = ind2.gameticFreq()
        self.assertDictEqual(gam_freq,{'AB': 0.24, 'Ab': 0.16000000000000003, 'aB': 0.36, 'ab': 0.24})
        
    def test_haploid(self):
        ind3 = pop.Individual('nombre_ind','nombre_pob',1,1,0,freq={'A':(0.4,0.6),'B':(0.6,0.4)},
                          d=0,R=0.5,mu=(0.1,0.1))
        #chromosomes 
        chr = ind3.chromosome
        self.assertEqual(list(chr.keys()),['c1'])
        genotype = ind3.genotype
        self.assertEqual(len(genotype['A']),1)
if __name__ == '__main__':
    unittest.main()
    