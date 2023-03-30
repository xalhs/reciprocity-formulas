import sympy
from sympy import *
import sys
from math import gcd
import itertools
from sympy.physics.quantum import TensorProduct
import copy
from src.util import *
from src.cokernel import compute_cokernel

def reciprocity_formula(mat , K_mat, tracesum = False):
    n = len(mat.row(0)) #dimension of matrix (nxn)
    m = len(K_mat.row(0)) #dimension of matrix (mxm)
    data_mat = compute_cokernel(mat)
    data_k_mat = compute_cokernel(K_mat)
    Q = data_mat[2] #the (reduced) linking form
    Q_K = data_k_mat[2]
    p_list = data_mat[0]  # Get the components of thetorsion homology
    k_list = data_k_mat[0]
    mat_reduced_det = data_mat[4] #determinant of the invertible part of the matrix
    K_mat_reduced_det = data_k_mat[4]
    l = len(p_list)  # dimension of torsion homology
    r = len(k_list)
    p_ran_tp = tuple([range(p_list[i]) for i in range(len(p_list))])
    k_ran_tp = tuple([range(k_list[i]) for i in range(len(k_list))])
    res1 = 0
    res2 = 0
    sgn1 = signature(mat) #signature of the matrix
    sgn2 = signature(K_mat)
    if tracesum:  #adds with the trace method if selected
        for vec1_tp in itertools.product(*k_ran_tp, repeat = n):
            vec1_new = Matrix(n , r , list(vec1_tp))
            res1 +=abs(K_mat_reduced_det)**(-l/2)*exp(pi*I*Trace(vec1_new*Q_k*vec1_new.T*mat).simplify())
        for vec2_tp in itertools.product(*p_ran_tp, repeat = m):
            vec2_new = Matrix(m, l , list(vec2_tp))
            res2 += abs(mat_reduced_det)**(-r/2)*exp(pi*I*sgn1*sgn2/4)*exp(-pi*I*Trace(vec2_new*Q*vec2_new.T*k_mat).simplify())
    else:     #adds with the normal method (default)
        for vec1_tp in itertools.product(*k_ran_tp, repeat = n):
            vec1 = Matrix(list(vec1_tp))
            res1 += abs(K_mat_reduced_det)**(-(n-mat.rank()/2))*exp(pi*I*(vec1.T*TensorProduct(mat, Q_K)*vec1)[0])

        for vec2_tp in itertools.product(*p_ran_tp, repeat = m):
            vec2 = Matrix(list(vec2_tp))
            res2 += abs(mat_reduced_det)**(-(m-K_mat.rank()/2))*exp(pi*I*sgn1*sgn2/4)*exp(-pi*I*(vec2.T*TensorProduct(K_mat, Q)*vec2)[0])

    print(res1.evalf())
    print(res2.evalf())

def rhs_sum(mat , K_mat):  #same as above
    n = len(mat.row(0))
    m = len(K_mat.row(0))
    data_mat = compute_cokernel(mat)
    data_k_mat = compute_cokernel(K_mat)
    Q = data_mat[2]
    Q_K = data_k_mat[2]
    p_list = data_mat[0]
    k_list = data_k_mat[0]
    mat_reduced_det = data_mat[4]
    K_mat_reduced_det = data_k_mat[4]
    l = len(p_list)
    r = len(k_list)
    p_ran_tp = tuple([range(p_list[i]) for i in range(len(p_list))])
    k_ran_tp = tuple([range(k_list[i]) for i in range(len(k_list))])
    res2 = 0
    sgn1 = signature(mat)
    sgn2 = signature(K_mat)
    for vec2_tp in itertools.product(*p_ran_tp, repeat = m):
        vec2 = Matrix(list(vec2_tp))
        res2 += abs(mat_reduced_det)**(-(m-K_mat.rank()/2))*exp(-pi*I*(vec2.T*TensorProduct(K_mat, Q)*vec2)[0])

    print(res2.evalf())
