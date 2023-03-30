import sympy
from sympy import *
from math import gcd

def rang(n , m = 0 ): #Redefining range to work better with negative numbers
    if m !=0:
        return range(n,m)
    if n >=0:
        return range(n )
    if n< 0:
        return range(0, n , -1)

def sub_matrix(mat , n):
    new_mat = Matrix(mat.col(0)[:n])
    for i in range(1, n):
        new_mat = new_mat.col_insert(i , Matrix(mat.col(i)[:n]))
    return new_mat

def K2(K , n , m , plus = True): #The second Kirby move, that adds row m to row n
    if n == m:
        raise("whaterror " + str(n) + " " + str(m))
    if plus == False:
        mult = -1
    else:
        mult = 1
    K[n,:] += mult*K[m,:]
    K[:,n] += mult*K[:,m]

def signature(mat): #The signature of the matrix, number of positive minus number of negative eigenvals
    sgn = 0
    for key in mat.eigenvals():
        sgn += sign(key)*mat.eigenvals()[key]
    return sgn

def set_of_points_within(mat):  #The points that are within the unit cell after the matrix acts on it
    if det(mat) == 0:
        print("matrix needs to be invertible for this algorithm")
        return False
    set_of_points = []
    n = len(mat.row(0))
    basis = []
    inv_basis = []
    for i in range(n):
        basis.append(eye(n).col(i))
        inv_basis.append(mat**-1*eye(n).col(i))

    searching = True
    coef = zeros( n , 1)
    set_of_points.append(mat*mat**-1*coef)
    while searching:
        for i in range(n):
            if not mat*vector_mod_Z(mat**-1*coef+inv_basis[i]) in set_of_points:
                coef[i] += 1
                set_of_points.append(mat*vector_mod_Z(mat**-1*coef))
                break
            else:
                coef[i] = 0
        if coef == zeros(n,1):
            if len(set_of_points) != abs(det(mat)):
                print("#"*50 + " Problem")
            return set_of_points

def split_matrix(mat):   #Separates invertible from non-invertible part of the matrix
    n = len(mat.row(0))
    basis = []
    for i in range(n):
        basis.append(eye(n).col(i))
    separating = True
    while separating:
        separating = False
        for nullvec in mat.nullspace():
            if nullvec in basis:
                continue
            else:
                separating = True
                denom = 1
                max_el = 0
                for comp in nullvec:
                    denom = lcm(denom , comp.q)
                    if abs(comp) > max_el:
                        max_el = abs(comp)
                nullvec = denom*nullvec
                if denom == 1:
                    glob_com_div = max_el
                    for comp in nullvec:
                        glob_com_div = gcd(glob_com_div , comp)
                    nullvec = nullvec / glob_com_div
                last_comp = 0
                init_index = n
                for i , comp in reversed(list(enumerate(nullvec))):
                    if last_comp ==0:
                        last_comp = -comp
                        init_index = i
                    else:
                        if comp == 0:
                            continue
                        while comp%last_comp != 0:  #n = last_comp
                            if last_comp == comp:  #m = comp
                                K2(mat , init_index , i , False)
                                comp -= last_comp
                            elif last_comp == -comp:
                                K2(mat , init_index , i , True)
                                comp += last_comp
                            elif abs(comp) > abs(last_comp):
                                if comp*last_comp>0:
                                    comp -= last_comp
                                    K2(mat , init_index , i , False)
                                else:
                                    comp += last_comp
                                    K2(mat , init_index , i)
                            elif abs(comp) < abs(last_comp):
                                if comp*last_comp>0:
                                    last_comp -= comp
                                    K2(mat , i , init_index , False)
                                else:
                                    last_comp += comp
                                    K2(mat , i , init_index)
                        ratio = abs(int(comp / last_comp))
                        for j in range(ratio):
                            if comp*last_comp > 0:
                                comp -= last_comp
                                K2(mat , init_index , i , False)
                            else:
                                comp += last_comp
                                K2(mat , init_index , i)
            break
    dim_nul = len(mat.nullspace())
    for i in range(dim_nul):
        if basis[n-i-1] not in mat.nullspace():
            for j in reversed(range(n-i-1)):
                if basis[j] in mat.nullspace():
                    mat.row_swap(j , n-i-1)
                    mat.col_swap(j, n-i-1)
    return mat

def is_mat_int(mat): #checks if all the components of a matrix are integers
    for i in mat:
        if not sympify(i).is_integer:
            return False
    return True

def vector_magnitude(mat): #gives the Euclidean norm of a vector in matrix form
    mag = 0
    for i in mat:
        mag += i**2
    return mag

def vector_mod_Z(mat): #gives the fractional part of the components of the vector in matrix form
    for i,comp in enumerate(mat):
        mat[i] = comp - floor(comp)
    return mat

def lcm(a,b): #least common multiple function because it does not exist in math package
    return int(a*b/gcd(a,b))
