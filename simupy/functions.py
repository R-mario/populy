
def outer2_product(a, b,iteration,finalDict):
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

def rename(i,j,iteration,k,finalDict):
    '''Set the proper letters according to the outer product'''
    # check if it is the first round
    if iteration==1:
        letrai =chr(ord('A')+iteration-1)
        #check if its the first position
        if i==1:
            letrai = letrai.lower()
    else:
        letrai = list(finalDict.keys())[k]

    if j==1:
        letraj = chr(ord('A')+iteration).lower()
    else:
        letraj = chr(ord('A')+iteration)
    return letrai,letraj
