
import pandas as pd 
import matplotlib.pyplot as plt 

from population import Population

class Superpop():
    '''
    This class allows to perform operations on serveral populations at once
    '''
    
    def __init__(self,kwargs):
        pass
    
    def createPops(popsize,n):
        """_summary_

        Args:
            popsize (_type_): _description_
            n (_type_): _description_

        Returns:
            _type_: _description_
        """
        # creamos
        create_sp = [Population(popsize)for x in range(n)]
        
        return create_sp
        
    def initPops(superP: list[Population]):
        printInfo = False
        # inicializamos
        [x.generateIndividuals()  for x in superP]
        #evolucionamos
        [x.evolvePop(printInfo=False) for x in superP]
        
        
        freqs = [d.alleleFreqs for d in superP]
        
        return freqs
    
    def plotPops(freqs):
        print(freqs)
        df = pd.DataFrame(freqs[0])
        fig,ax=plt.subplots(1,2)
        popListnames = ['p'+str(i) for i in range(len(freqs))]
        for i in range(len(freqs)):
            ax[0].plot(freqs[i]['A'])
            ax[1].plot(freqs[i]['B'])

        plt.show()
        
        
        
        
        
        
if __name__ == '__main__':
    
    superpop = Superpop.createPops(100,n=20)
    
    alFreq_sp = Superpop.initPops(superpop)
    
    
    Superpop.plotPops(alFreq_sp)
        