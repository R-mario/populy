import unittest
from simupy.individual import Individual


class individualTest(unittest.TestCase):

    def test_individual_activation(self):
        global ind1
        # instanciamos

        ind1 = Individual('a','test',1,0,0,freq={'A':(0.4,0.6),'B':(0.6,0.4)},
        d=0,R=0.5,mu=(0.1,0.1))
        # comprobamos que devuelve un string
        self.assertTrue(ind1)
    def gameticTest(self):
        chromosome = ind1.gameticFreq()
        self.assertDictEqual(chromosome,{'AB': 0.24, 'Ab': 0.16000000000000003, 'aB': 0.36, 'ab': 0.24})

if __name__ == '__main__':
    unittest.main()
    