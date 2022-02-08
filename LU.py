import numpy as np

#####2######
#a#
def LU_Dec(matrix):
    N = np.shape(matrix)[0] #number of rows(or columns)
    U = np.zeros((N,N))
    L = np.zeros((N,N))
    
    for j in range(N):
         
        L[j][j] == 1 #fixes the diagonal values to 1
        
        for i in range(j+1): #since j is from 0 to N-1, need to add 1 for range to N
            sum1= sum(L[i][k]*U[k][j] for k in range(i))
            U[i][j]=matrix[i][j]-sum1
        
            
        for i in range(j, N): #i is from j to N 
            sum2= sum(L[i][k]*U[k][j] for k in range(j))
            L[i][j]=(1/U[j][j])*(matrix[i][j]-sum2)
            
    print('Upper is', U)        
    print('Lower is', L)
    #print('Product is upper times lower:', U,'X',L)
    return L,U
#b#                
Dec=LU_Dec([[3,1,0,0,0],[3,9,4,0,0],[0,9,20,10,0],[0,0,-22,31,-25],[0,0,0,-55,60]])
'''LU decomposition is 
  U=
 [  3.           1.           0.           0.           0.        ]
 [  0.           8.           4.           0.           0.        ]
 [  0.           0.          15.5         10.           0.        ]
 [  0.           0.           0.          45.19354839 -25.        ]
 [  0.           0.           0.           0.          29.57530335]
  L=
 [ 1.          0.          0.          0.          0.        ]
 [ 1.          1.          0.          0.          0.        ]
 [ 0.          1.125       1.          0.          0.        ]
 [ 0.          0.         -1.41935484  1.          0.        ]
 [ 0.          0.          0.         -1.21698787  1.        ]
'''
def det(decomposed_matrix):
    N = np.shape(decomposed_matrix[1])[0] #Second element of decomposed matrix is 
                                          #upper matrix, from which we need the 
                                          #number of rows, which is first element
                                          #of the returned value from shape function
    
    V = np.zeros((N,N)) #creating a list to be used later to multiply
                                  #diagonal elements in U
    det = 1
    for i in range(N):
        V[i][i] = decomposed_matrix[1][i][i]
        det *= V[i][i]
    
    print('Determinant is', det)
    return det

det(Dec) 

'''Determinant is 497220.0'''

#c#
def matrix_eqn_solver(L,U,b): #LUx=b (Vector eqn)
    N = np.shape(b)[0]
    
    x, y = np.zeros((N,N)), np.zeros((N,N))
    
    for i in range(N):
        #forward substitution
        sum1 = sum(L[i][j]*y[j] for j in range(i))
        y[i]=(1/L[i][i])*(b[i]-sum1)
            
            
        
    for i in range(N-1,-1,-1):
        #backward substitution
        sum2 = sum(U[i][j]*x[j] for j in range(i, N))
        x[i] = (1/U[i][i])*(y[i]-sum2)
            
        
    #print ('The solution is', 'x=', x)
    
    return x

#d#
matrix_eqn_solver(Dec[0], Dec[1], [2,5,-4,8,9])
'''x= [0.4561763404529182, 0.6314709786412454, -0.5129419572824906, 0.05756003378786046, 0.20276336430553876]'''

#e#
#we want to find A^-1 for AA^-1=I so we set A^-1 as x and I as b
inverse=matrix_eqn_solver(Dec[0], Dec[1], np.identity(5))
'''Inverse of A is 
[[ 0.37938941 -0.04605607  0.00390169 -0.00482684 -0.00201118]
 [-0.13816822  0.13816822 -0.01170508  0.01448051  0.00603355]
 [ 0.02633643 -0.02633643  0.02341016 -0.02896102 -0.01206709]
 [ 0.07167853 -0.07167853  0.06371425  0.04488959  0.01870399]
 [ 0.06570532 -0.06570532  0.05840473  0.04114879  0.03381199]]
'''
