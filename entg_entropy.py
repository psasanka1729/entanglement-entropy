#!/usr/bin/env python
# coding: utf-8

import numpy as np
import random


# Number of qubits.
N = 4

L = N // 2 # Length of half cut number of qubits.

'''
    Intiating a wave function with a lsit of size 2**N with all element as zeros.

''' 

Psi_List = [0]*(2**N)




'''
    It is not necessary to input a normalized wave function. 
    It will be normalized later so that trace(rho) = 1.

'''


'''
    Enter the non-zero coefficients of the wavefunction psi.

''' 

Psi_List[0] = 1/2
Psi_List[3] = -1/2
Psi_List[12] = 1/2
Psi_List[15] = -1/2






'''
    The following function takes a wavefunction as input and returns its entropy.

'''

def Entropy(Wavefunction):




    # Converting the list to a numpy matrix.
    Psi = np.matrix(Wavefunction).reshape(len(Wavefunction),1) # Psi column matrix.

    # Normalizing Psi.
    Psi = Psi/np.linalg.norm(Psi)


      
    
    def psi(s):
        return Psi[(2**L)*s:(2**L)*s + 2**L]   
    
      
    '''
        psi(s_p) is a row matrix/vector. psi(s) is a column matrix/vector.      
        Dimension of rhoA is N/2 x N/2. 
        The element <s|rhoA|sp> is given by psi_sp^\dagger * psi_s.
        
    ''' 

    def rhoA(s,s_p): # <s|rho_A|s_p>

        # psi(s_p)^\dagger * psi(s) is the element of (s,s_p) of rho_AB.  
        return psi(s_p).getH() * psi(s)
    
    
    
    def rhoA_Matrix(N):
    
        M = np.zeros((N,N)) # 0 to N-1.
    
        '''
            rho is Hermitian, it is sufficient to calculate the elements above the diagonal.
            The the elements below the diagonal can be replace by the complex cpnjugate of the
            elements above the diagonal.
        '''
        for i in range(N):
            for j in range(N):
            
                if i <= j : # Above the diagonal (i,j) i<j.
                
                    M[i,j] = rhoA(i,j)[0,0]
                
                else: # Below the diagonal (i,j) i>j.
                
                    M[i,j] = np.conjugate(M[j,i])
        return M    
    
    
    '''
        w is the diagonal of the diagonalized matrix rhoA.

    '''
    
    w, v = np.linalg.eig(rhoA_Matrix(N))
    
    w = w.real

    '''
        The following loop calculates S = - sum \lamba_i * log(\lambda_i).

    '''
    
    DL = np.zeros(N) # Creating an array for log w with zeros.
    
    for i in range(len(w)):
    
        if abs(w[i]) < 1.e-8: # log of zero gives nan.
        
            pass # Leave the log(zero) element as zero.
    
        else:
        
            DL[i] = np.log(w[i])
        
    # Entropy = -Tr(rho * log(rho)).        
    return -sum(w*DL)






def Bin2Dec(BinaryNumber): # Converts binary to decimal numbers.
    return int(str(BinaryNumber),2)


def Dec2Bin(DecimalNumber): # Converts decimal to binary numbers.
    return bin(DecimalNumber).replace("0b", "")



List = [i for i in range(2**N)] 


'''
The following function converts all numbers in decimals in the above list  from 0 to 2^N -1 to binary.

''' 
def List_Bin(List):
    
    l = []
    
    for i in List:
        
        i_Bin = Dec2Bin(i)
              
        
        '''
        While converting numbers from decimal to binary, for example, 1 is mapped to 1, to make sure that
        every numbers have N qubits in them, the following loop adds leading zeros to make the
        length of the binary string equal to N. Now, 1 is mapped to 000.....1 (string of length N).
        
        '''
        
        while len(i_Bin) < N: 
            
            i_Bin = '0'+i_Bin # This loop adds leading zeros.
            
        l.append(i_Bin)
        
    return l





'''
    The following function takes a binary string as input and rolls the qubits by one and
    returns the rolled string.

'''

def Roll_String(Binary_String):
    
    return Binary_String[-1] + Binary_String[:-1]







'''
    The following function takes a wavefunction as input and performs one roll on the qubits and
    returns the resultant wavefunction.

'''

def Psi_Roll(Inital_Psi):
    
    
    
    '''
        The following list contains all possible 2^N qubits after one roll is performed on them.
        For example, the first position 0001 is changed to 1000.
    
    '''
    
    Rl = [Roll_String(i) for i in List_Bin(List)] # Rolls every string in the list List by one qubit.

   

    
    ''' 
        The following list contains the previous list but in decimal numbers. For example,
        for N =4, the first position 1 is changed to 8.
        
    
    '''
    
    Rl_d = [Bin2Dec(i) for i in Rl] # Converts the rolled binary string to decimal number.


    '''
        The following loop rearranges the coefficients of Psi after rolling. For example, for N = 4,
        if the first coefficient 0001 is mapped to the eighth coefficient 1000 after one rotation of
        the qubits. The coefficient of the rolled Psi in the i ^ th position is in the Rl_d[i] ^ th positon
        of the initial Psi.
    
    '''
    
    
    Psi_Rolled = []

    for i in range(2**N): 
    
        Psi_Rolled.append(Inital_Psi[Rl_d[i]]) # Rearranging the coefficients according to the list l_d.
        
    return Psi_Rolled






'''
    The following function performs specified number of rolls Num on the qubits.

'''

def N_Rolled(Num, Initial_Psi): # Use loop only for postive N.
    
    s = Psi_Roll(Initial_Psi) # One roll.
    
    for i in range(Num-1): # Loop performing remaining N-1 rolls.
        
        s = Psi_Roll(s)
        
    return np.matrix(s).reshape(2**N,1) # Returning the rolled wavefunction as a matrix.



print('Entropy after 0 roll = ', Entropy(Psi_List)) # Zero roll case.

for i in range(1,N): # Prints rolls from 1 to N-1.
    print('Entropy after '+ str(i) + ' roll =', Entropy(N_Rolled(i,Psi_List)))

