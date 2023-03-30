import sympy
from sympy import *
import sys
from math import gcd
import itertools
from sympy.physics.quantum import TensorProduct
import copy
from src.util import *

def phase_calc(res): #approximate function to get the phase of a number
    if abs(re(res)) < 0.001:   #works best if the phase is quantized like in gauss sums
        ang = pi/2
        if abs(im(res)) < 0.001:
            ang = 0
        elif im(res) < 0:
            ang = -ang
    elif re(res) <0:
        if abs(im(res)) < 0.001:
            ang = pi
        elif im(res) < 0:
            ang = -(pi - abs(atan(im(res)/re(res))))
        else:
            ang = (pi - abs(atan(im(res)/re(res))))
    elif re(res) > 0:
        if abs(im(res)) < 0.001:
            ang = 0
        elif im(res) < 0:
            ang = -abs(atan(im(res)/re(res)))
        else:
            ang = abs(atan(im(res)/re(res)))
    return ang

def phase(res): #another way to get the phase
    return sympy.log( res ).as_real_imag()[1]

def fraction_of_pi(a): #gives the closest multiple of pi/4 of a number that is supposed to be n*pi/4
    if any(ext in str(a) for ext in ["atan" , "log"]):
        for o in rang(-3,5):
            if abs(a-o*pi/4).evalf()<0.001:
                return o*pi/4
    return a

def ordering(list): #rewrites a list of integers so that they are multiples of each other while keeping their product the same.
    for i in range(len(list)):
        for j in range(i+1 , len(list)):
            g_c_d = gcd(list[i] , list[j])
            l_c_m = lcm(list[i] , list[j])
            list[i] = g_c_d
            list[j] = l_c_m
    return list

def pow_of_2_in_number(num): #gives the greatest power of 2 in a positive integer
    if num >0:
        power = 2
        while True:
            if num%power != 0:
                return int(power/2)
            else:
                power *= 2
    else:
        print("It was made to calculate positive integers")

def Jacobi_sym(q,p): #Jacobi symbol function after eliminating powers of 2 in p
    p = int(p/pow_of_2_in_number(p))
    return jacobi_symbol(q,p)

def gauss_sum(p,q): #gives the regular gauss sum exp(i*pi*q*n^2/p) from n=0 to p-1
    res = 0
    for n in rang(p):
        res += exp(q*pi*I*n*n/(p))
    return res

def is_within(matrix , point): #checks if a point is within the parallelepiped created by the matrix
    for el in matrix**-1*point.T:
        if not (0 <= el and el < 1):
            return False
    return True

def orth_diag(mat): #algorithm for orthogonal diagonalization
    n = len(mat.row(0))
    iden = eye(n)
    for k in range(n):
        if mat[k,k]== 0:
            for i in rang(k+1 , n):
                if mat[k,i]!=0:
                    mat[k,:]+= mat[i,:]
                    iden[k,:] += iden[i,:]
                    mat[:,k] += mat[:,i]
        for i in range(k+1,n):
            m1 = mat[k,k]
            m2 = mat[i,k]
            mat[i,:] = -m1*mat[i,:]+m2*mat[k,:]
            iden[i,:] = -m1*iden[i,:]+m2*iden[k,:]
            mat[:,i] = -m1*mat[:,i]+m2*mat[:,k]
            print([mat,iden])

    return [mat,iden]

def linking_mat_of_lens_space(p,q): #gives an integer linking matrix of lens space L(p,q) (saveliev)
    frac_exp = list(continued_fraction_iterator(Rational(p,q)))
    n = len(frac_exp)
    mat = zeros(n)
    for i in range(n):
        for j in range(n):
            if i==j:
                mat[i,i] = -frac_exp[i]*(-1)**i
            elif abs(i-j) == 1:
                mat[i,j] = 1
    return mat

def tors_of_cont_frac(p,q): #gives the torsion numbers for the linking matrix of L(p,q) in form [p,1,...,1]
    frac_exp = list(continued_fraction_iterator(Rational(p,q)))
    n = len(frac_exp)
    lst = []
    for i in range(n):
        lst.append(1)
    lst[0] = p
    return lst

def reciprocity1_extra(mat , p_list , k): #two reciprocity formulas, they also take as input the list of torsion coefficients
    p_ran_tp = tuple([rang(p_list[i]) for i in rang(len(p_list))])
    l = len(mat.row(0))
    sg = signature(mat)
    res1 = 0.0
    res2 = 0.0
    for vec1_tp in itertools.product(rang(2*k) , repeat = l):
        vec1 = Matrix(list(vec1_tp))
        res1 += (2*k)**(-l/2)*exp(pi*I*(1/(2*k))*(vec1.T*mat*vec1)[0])
    for vec2_tp in itertools.product(*p_ran_tp):
        vec2 = Matrix(list(vec2_tp))
        res2 += exp(pi*I*sg/4)*abs(det(mat))**(-1/2)*exp(-pi*I*2*k*(vec2.T*mat**-1*vec2)[0])

    print(res1.evalf())
    try:
        print(res2.evalf())
    except:
        print(res2)

def reciprocity2_extra(mat , k_mat , p_list , k_list): #the first formula is for 1x1 dim K while the second for arbitrary
    p_ran_tp = tuple([rang(p_list[i]) for i in rang(len(p_list))])
    k_ran_tp = tuple([rang(k_list[i]) for i in rang(len(k_list))])
    l = len(mat.row(0))
    r = len(k_mat.row(0))
    res1 = 0
    res2 = 0
    sgn1 = signature(mat)
    sgn2 = signature(k_mat)
    for vec1_tp in itertools.product(*k_ran_tp, repeat = l):
        vec1 = Matrix(list(vec1_tp))
        res1 += det(k_mat)**(-l/2)*exp(pi*I*(vec1.T*TensorProduct(mat, (k_mat)**-1)*vec1)[0])
    for vec2_tp in itertools.product(*p_ran_tp , repeat = r):
        vec2 = Matrix(list(vec2_tp))
        res2 += (det(mat))**(-r/2)*exp(pi*I*sgn1*sgn2/4)*exp(-pi*I*(vec2.T*TensorProduct(k_mat, mat**-1)*vec2)[0])
    print(res1.evalf())
    print(res2.evalf())
