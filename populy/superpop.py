
import pandas as pd 
import matplotlib.pyplot as plt 
import os
import numpy as np
from population import Population

class Superpop():
    '''
    This class allows to perform operations on serveral populations at once
    '''
    
    def __init__(self,popsize,n,**kwargs):
        '''
        Args:
            popsize (int): tamaño de cada población
            n (int): número de poblaciones
            kwargs: los mismos de poblacion
            '''

        self.popsize = popsize
        self.n = n
        self.sPop = [Population(popsize,**kwargs)for x in range(n)] 
    
        
    def startPops(self,gens,**kwargs):
        cada = 5
        for key,value in kwargs.items():
            if key == 'every':
                cada = value

            
        printInfo = False
        # inicializamos
        [x.initIndividuals()  for x in self.sPop]
        #evolucionamos
        [x.evolvePop(gens,every=cada,printInfo=False) for x in self.sPop]
        # guarda las frecuencias alelicas acumuladas en una lista
        self.freqs = [d.f_ale_acc for d in self.sPop]
        os.system("cls")
        return self.freqs
    
    def plotPops(self):
    
        # numero de generaciones
        gens = self.sPop[0].gen
        # recogida de informacion
        steps = self.sPop[0].steps
        
        labels = ['gen.'+str(x) for x in range(0,gens+1,steps)]
        
              
        popListnames = ['p'+str(i) for i in range(len(self.freqs))]
        
        fig,ax=plt.subplots(1,2,figsize=(13,6))

        for j,let in enumerate(self.freqs[0].keys()):
            for i in range(len(self.freqs)): 
                ax[j].set_title(f'frecuencia alelica de {let}')
                ax[j].plot(self.freqs[i][let])
                
        new_steps = int(gens/20) if len(labels)>10 else 1
        plt.setp(ax, xticks=range(0,len(labels),new_steps), xticklabels=labels[::new_steps])
        plt.show()
        
        
if __name__ == '__main__':
    
    superpop = Superpop(100,n=5,ploidy=2,fit={'A':(0.9,1)},rnd=False)
    
    alFreq_sp = superpop.startPops(gens=100)
    
    
    superpop.plotPops()
        