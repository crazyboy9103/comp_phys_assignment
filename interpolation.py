import numpy as np
import numpy.fft as nft
import matplotlib.pyplot as plt
from LU import LU_Dec
from LU import matrix_eqn_solver
#############3########################
#a)

#use of linear interpolation
#x represents independent variable, y dependent variable
def linear_inter(x,y):
    x_inter = []
    y_inter = []
    
    for i in range(len(y)):
        if i+1 != len(y):
            z = np.linspace(x[i], x[i+1], 10)
            
            A = (x[i+1]-z)/(x[i+1]-x[i])
            B = 1-A
            f = A*y[i]+B*y[i+1]
            
            #print('The function evaluates as', f, 'at point x=', z)
            
            x_inter.extend(z)
            y_inter.extend(f)
        
        else:
            continue
            
    inter_function = [x_inter, y_inter]
                
    return inter_function

#b)
class cubic_spline():
    def __init__(self, x, y):
        self.__x = x
        self.__y = y
        
    def coefficient_ac_i(self,i):
        #The equation is a_i f''_(i-1)+b_i f''_(i)+c_i f''_(i+1)=F_i
        a_i = (self.__x[i]-self.__x[i-1])/6
        c_i = (self.__x[i+1]-self.__x[i])/6
        
        return a_i, c_i
    
    def coefficient_bF_i(self, i):
        b_i = (self.__x[i+1]-self.__x[i-1])/3
        F_i = ((self.__y[i+1]-self.__y[i])/(self.__x[i+1]-self.__x[i]))-((self.__y[i]-self.__y[i-1])/(self.__x[i]-self.__x[i-1]))
        
        return b_i, F_i
    
    def matrix_eqn_setup(self):
        N = np.shape(self.__x)[0]
        #for N points, there are N-1 intervals, or equivalently, coefficients
        A = np.zeros((N-1,N-1))
        F = np.zeros((N-1,N-1))
        
        for i in range(N-2):
            A[i][i+1] = self.coefficient_ac_i(i)[1]
            A[i+1][i] = self.coefficient_ac_i(i+1)[0]
        
        for i in range(N-1):
            A[i][i]=self.coefficient_bF_i(i)[0] 
            F[i][0]=self.coefficient_bF_i(i)[1] #set F as a vector inside a matrix
            
        return A, F
    
    def matrix_solution(self):
        #print('A is', self.matrix_eqn_setup()[0])
        #print('F is', self.matrix_eqn_setup()[1])
        dec_matrix = LU_Dec(self.matrix_eqn_setup()[0])
        
        #print('L is', dec_matrix[0])
        solution = matrix_eqn_solver(dec_matrix[0], dec_matrix[1], self.matrix_eqn_setup()[1])
        
        n = len(solution)
        #print('The solution is', solution)
        solution_vector = np.zeros(n)
    
        for i in range(n):
            solution_vector[i] += solution[i][0]
        
        #print('The solution vector is', solution_vector)
        return solution_vector  
                    
    def computing_function(self):
        x = []
        y = []
        for i in range(len(self.__y)-1):
            if i+1 != len(self.__y)-1:
                z = np.linspace(self.__x[i], self.__x[i+1],10)
                A = (self.__x[i+1]-z)/(self.__x[i+1]-self.__x[i])
                B = 1-A
                C = 1/6 * (A**3-A)*((self.__x[i+1]-self.__x[i])**2)
                D = 1/6 * (B**3-B)*((self.__x[i+1]-self.__x[i])**2)
                
                second_der = self.matrix_solution
                f = A*self.__y[i]+B*self.__y[i+1]+C*second_der()[i]+D*second_der()[i+1]
                x.extend(z)
                y.extend(f)
            
            else:
                continue
        
        function = [x,y]
        
        return function
   
#%%
#c#
table = [[-2.1, -1.45, -1.3,-0.2, 0.1, 0.15, 0.8, 1.1, 1.5, 2.8, 3.8],\
[0.012155, 0.122151, 0.184520, 0.960789, 0.990050, 0.977751, 0.527292, 0.298197, 0.105399, 3.936690E-4, 5.355348E-7]]

linear = linear_inter(table[0], table[1])
cubic_interpolated = cubic_spline(table[0], table[1]).computing_function()
plt.figure(3)
plt.plot(linear[0], linear[1],'-r')
plt.plot(cubic_interpolated[0], cubic_interpolated[1], '-b')
