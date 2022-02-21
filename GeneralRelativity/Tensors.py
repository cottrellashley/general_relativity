import sympy as smp
import numpy as np
import itertools as it

       
class Tensors:
    
  
    def __init__(self, A, B, SummedIndices):
        self.A = A
        self.B = B
        self.SummedIndices = SummedIndices
        self.A_dimention = int(smp.shape(A)[0])
        self.B_dimention = int(smp.shape(B)[0])
        self.A_number_of_indices = len(smp.shape(A))
        self.B_number_of_indices = len(smp.shape(B))
        self.A_summed_indices = list(it.product(np.arange(0, self.A_dimention, 1), repeat = self.A_number_of_indices))
        self.B_summed_indices = list(it.product(np.arange(0, self.B_dimention, 1), repeat = self.B_number_of_indices))
        self.Answer_number_of_indices = (self.A_number_of_indices + self.B_number_of_indices) - 2*len(self.SummedIndices)
 
 
    def Shape_of_Answer(self):
        D = self.A_dimention
        N = self.Answer_number_of_indices
        Ans_shape = ()
        for i in range(N):
            y = list(Ans_shape)
            y.append(D)
            Ans_shape = tuple(y)
        return Ans_shape

       
    def Answer_Summed_Indices(self):
        D = self.A_dimention
        Shape = self.Shape_of_Answer()
        return list(it.product(np.arange(0, D, 1), repeat = len(Shape)))
 
 
    def Answer_Tensor(self):
        D = self.A_dimention
        A = int(self.Answer_number_of_indices)
        N = int(self.Answer_number_of_indices)
        Shape = self.Shape_of_Answer()
        return smp.MutableDenseNDimArray(smp.zeros(D**A), Shape)

       
    def Left(self):
        I = self.SummedIndices
        A = len(smp.shape(self.A))
        B = len(smp.shape(self.B))
        a1 = np.arange(0, A, 1)
        b1 = np.arange(0, B, 1)
        Ia = []
        Ib = []
        for i in range(len(I)):
            Ia.append(I[i][0])
            Ib.append(I[i][1])
        return np.delete(a1,Ia).tolist(), np.delete(b1,Ib).tolist()

       
    def Input(self):
        n = self.Left()
        I1f, I2f = [ [] , [] ]
        for j in range(len(n[0])):
            I1f.append([j, n[0][j]])
        for i in range(len(n[0]), len(n[1]) + len(n[0]) ):
            I2f.append([i, n[1][i-len(n[0])]])
        return [I1f, I2f]

    
    def bitwise_and(self,i,j):
        I = self.SummedIndices
        a = self.A_summed_indices
        b = self.B_summed_indices
        c = self.Answer_Summed_Indices()
        B = [a[i][I[0][0]] == b[j][I[0][1]]]
        for k in range(len(I) - 1):
            B.append(B[k] and a[i][I[k+1][0]] == b[j][I[k+1][1]])
        return B[-1]


    def ans1_bitwise_and(self,i,j):
        I = self.Input()[1]
        a = self.A_summed_indices
        b = self.B_summed_indices
        c = self.Answer_Summed_Indices()
        if not len(I) == 0:
            B = [c[i][I[0][0]] == b[j][I[0][1]]]
            for k in range(len(I) - 1):
                B.append(B[k] and c[i][I[k+1][0]] == b[j][I[k+1][1]])
            return B[-1]
        else:
            return print("Function is irrelevant for the final result.")


    def ans2_bitwise_and(self,i,j):
        I = self.Input()[0]
        a = self.A_summed_indices
        b = self.B_summed_indices
        c = self.Answer_Summed_Indices()
        if not len(I) == 0:
            B = [c[i][I[0][0]] == a[j][I[0][1]]]
            for k in range(len(I) - 1):
                B.append(B[k] and c[i][I[k+1][0]] == a[j][I[k+1][1]])
            return B[-1]
        else:
            return print("Function is irrelevant for the final result.")


    def bitwise_and_func(self, n, j, i):
        I1 = self.Input()[0]
        I2 = self.Input()[1]
        
        if int(len(I2)) != 0 and int(len(I1)) != 0:
            return self.ans1_bitwise_and(n,j) and self.ans2_bitwise_and(n,i)

        if int(len(I2)) != 0 and int(len(I1)) == 0:
            return self.ans1_bitwise_and(n,j)
                
        if int(len(I1)) != 0 and int(len(I2)) == 0:
            return self.ans2_bitwise_and(n,i)
    
    
    def Contraction(self):
        A = len(self.A_summed_indices)
        B = len(self.B_summed_indices)
        C1 = self.Answer_Summed_Indices()
        C = len(C1)
        A1 = self.A
        B1 = self.B
        C1 = self.Answer_Tensor()
        a = self.A_summed_indices
        b = self.B_summed_indices
        c = self.Answer_Summed_Indices()
        
        if len(C1) > 1:
            for n in range(C):
                for i in range(A):
                    for j in range(B):
                        if self.bitwise_and(i,j) and self.bitwise_and_func(n, j, i):
                            C1[c[n]] += B1[b[j]]*A1[a[i]]
            return C1
        
        if len(C1) == 1:
            Constant = float()
            for i in range(A):
                for j in range(B):
                    if self.bitwise_and(i,j):
                        Constant += B1[b[j]]*A1[a[i]]
            return Constant