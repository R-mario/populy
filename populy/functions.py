
import random

def outer_product(a: int, b : int,iteration :int,finalDict: dict):
    """Calculate the outer product of two vectors and stores them on a dictionary"""

    d = dict()

    values = list()
    names = list()
    k=0
    # go over first vector values
    for i,val_i in enumerate(a):
        #go over second vector values
        for j,val_j in enumerate(b):
                #append values into new list
                values.append(val_i*val_j)

                letrai,letraj = rename(i,j,iteration,k,finalDict)
                names.append(letrai+letraj)
                if j%2==1:
                    k+=1
    d = dict(zip(names,values))
    return d

def outer_product2(f):
    """Calculate the outer product of two vectors and stores them on a dictionary"""
    # lista de tuplas con los valores de frecuencia alelica
    fValues = list(f.values())
    # primera tupla
    a = fValues[0]
    finalD = dict()
    # recorremos todos los valores
    for x in range(1,len(fValues)):
        d = dict()
        values = list()
        names = list()
        k=0
        # recorremos la tupla numero i
        for i,val_i in enumerate(a):
            #recorremos el resto de tuplas
            for j,val_j in enumerate(fValues[x]):
                    # multiplicamos el valor i de la primera tupla
                    # por el valor j de la tupla x
                    values.append(val_i*val_j)
                    # cambiamos el nombre
                    letrai,letraj = rename(i,j,x,k,finalD)
                    names.append(letrai+letraj)
                    if j%2==1:
                        k+=1
        finalD = dict(zip(names,values))
        a = list(finalD.values())
    if len(list(f.keys()))==1: # si solo hay un locus
        # la letra (prob. A)
        letter = list(f.keys())
        # lista de keys (A,a)
        letter.append(letter[0].lower())
        # lista de valores ej: (0.4,0,6)
        values = list(f.values())[0]
        # se hace diccionario
        finalD = dict(zip(letter,values))
        
    return finalD

def rename(i,j,iteration,k,finalDict):
    '''Set the proper letters according to the outer product'''
    # check if it is the first round
    if iteration==1:
        letrai =chr(ord('A')+iteration-1)
        # check if its the first position
        if i==1:
            letrai = letrai.lower()
    else:
        letrai = list(finalDict.keys())[k]

    if j==1:
        letraj = chr(ord('A')+iteration).lower()
    else:
        letraj = chr(ord('A')+iteration)
    return letrai,letraj


def fitness(fit,child_gen):
    '''Fitness function
    0 = no selection,
    1 = dominant selection (1,0.9,0.91)
    2 = recessive selection (0.81,0.9,1)
    3 = double recessive selection 'aabb'=0
    '''
    def live_die(p):
        """_summary_

        Args:
            p (_type_): _description_

        Returns:
            _type_: _description_
        """
        randNum = random.random()
        if randNum < p:
            return True
        else:
            return False
  
    
    if fit==0:
        return True
    if fit==1:
        if 'Aa'== child_gen['A']:
            return live_die()
        if 'aa'== child_gen['A']:
            return live_die(0.81)
    if fit==2:
        if 'Aa'== child_gen['A']:
            return live_die(0.9)
        if 'AA'== child_gen['A']:
            return live_die(0.81)
    if fit==3:
        if 'aa'== child_gen['A'] and 'bb'== child_gen['B']:
            return False
    
    return True
        
        
    
