import sympy
from sympy import *
import sys
from math import gcd
import itertools
from sympy.physics.quantum import TensorProduct
import copy
from src.util import *

def classify_points(mat , list_of_points):
    dict_of_points = {}
    dict_of_inv_points = {}
    max_order = 0
    ord = 0
    max_ord_points = []
    for point in list_of_points:
        if tuple(point) in dict_of_points:
            continue
        [ord , new_dict] = order_of_point(mat ,point)
        for el in new_dict:
            dict_of_points[tuple(mat*Matrix(list(el)))] = new_dict[el]
            dict_of_inv_points[el] = new_dict[el]
        if ord >= max_order:
            new_max_points = [point]
            for el in new_dict:
                if new_dict[el] == ord:
                    new_max_points.append(mat*Matrix(list(el)))
        if ord == max_order:
            max_ord_points += new_max_points
        elif ord > max_order:
            max_order = ord
            max_ord_points = new_max_points

        dict_of_points[tuple(point)] = ord
        dict_of_inv_points[tuple(mat**-1*point)] = ord
    min_mag = 0
    for point in max_ord_points:
        vec_mag = vector_magnitude(point)
        if vec_mag ==0:
            min_vec = point
            break
        elif min_mag ==0 or vec_mag < min_mag:
            min_mag = vec_mag
            min_vec = point

    return [dict_of_points, max_order , max_ord_points , min_vec , dict_of_inv_points , mat**-1*min_vec]

def classify_new_classes(dict_of_classes_org, repr , new_equiv , order_of_equiv):
    dict_of_reprs = {}
    max_order = 0
    ord = 0
    max_ord_reprs = []
    for point in repr:
        if tuple(point) in dict_of_reprs:
            continue
        [ord , new_dict] = order_of_new_equiv_point(point , dict_of_classes_org,  new_equiv , order_of_equiv  )
        for el in new_dict:
            dict_of_reprs[el] = new_dict[el]
        if ord >= max_order:
            new_max_reprs = [point]
            for el in new_dict:
                if new_dict[el] == ord:
                    new_max_reprs.append(Matrix(list(el)))
        if ord == max_order:
            max_ord_reprs+= new_max_reprs
        elif ord > max_order:
            max_order = ord
            max_ord_reprs = new_max_reprs

        dict_of_reprs[tuple(point)] = ord
    min_mag = 0
    for point in max_ord_reprs:
        vec_mag = vector_magnitude(point)
        if vec_mag ==0:
            min_vec = point
            break
        elif min_mag ==0 or vec_mag < min_mag:
            min_mag = vec_mag
            min_vec = point
    return [dict_of_reprs , max_order , max_ord_reprs , min_vec]

def quotient_space(list_of_inv_points_org , list_of_equiv ,orders_of_equiv):
    list_of_inv_points = copy.deepcopy(list_of_inv_points_org)
    classes = []
    dict_classes = {}
    repr = []
    num_of_classes = 0
    eq_ran_tp = tuple([range(orders_of_equiv[i]) for i in range(len(orders_of_equiv))])
    while list_of_inv_points != []:
        classes.append([])
        starting_vec = list_of_inv_points[0]
        for comb in itertools.product(*eq_ran_tp):
            res_vec = starting_vec
            for i,mult in enumerate(list(comb)):
                res_vec += mult*list_of_equiv[i]
            res_vec = vector_mod_Z(res_vec)
            classes[num_of_classes].append(res_vec)
            list_of_inv_points.remove(res_vec)
            min_mag = 0
        for j in classes[num_of_classes]:
            vec_mag = vector_magnitude(j)
            if vec_mag ==0:
                min_vec = j
                break
            elif min_mag ==0 or vec_mag < min_mag:
                min_mag = vec_mag
                min_vec = j
        repr.append(min_vec)
        dict_classes[tuple(min_vec)] = classes[num_of_classes]
        num_of_classes += 1
    return [dict_classes , classes , repr ]

def new_quot_space(dict_of_classes_org ,repr_org , new_equiv ,order_of_equiv):
    dict_of_classes = copy.deepcopy(dict_of_classes_org)
    repr = copy.deepcopy(repr_org)
    new_classes = []
    dict_new_classes = {}
    classes_of_repr = []
    new_repr = []
    num_of_classes = 0
    while repr!= []:
        classes_of_repr.append([])
        new_classes.append([])
        point = repr[0]
        for i in range(order_of_equiv):
            new_point = vector_mod_Z(point + i*new_equiv)
            for repres in repr:
                for other_point in dict_of_classes[tuple(repres)]:
                    if new_point == other_point:
                        classes_of_repr[num_of_classes].append(repres)
                        rep_found = repres
            repr.remove(rep_found)
        min_mag = 0
        for j in classes_of_repr[num_of_classes]:
           #new_classes[num_of_classes].append(dict_of_classes[tuple(j)])
            new_classes[num_of_classes] += dict_of_classes[tuple(j)] #the above caused an infinite loop with larger matrices, this fixes the bug
            vec_mag = vector_magnitude(j)
            if vec_mag ==0:
                min_vec = j
                break
            elif min_mag ==0 or vec_mag < min_mag:
                min_mag = vec_mag
                min_vec = j
        new_repr.append(min_vec)
        dict_new_classes[tuple(min_vec)] = new_classes[num_of_classes]
        num_of_classes += 1
    return [dict_new_classes , new_classes , new_repr ]

def order_of_point(mat, point):
    order = 0
    other_points_dict = {}
    while True:
        order +=1
        if is_mat_int(mat**-1*(order*point)):
            for i in other_points_dict:
                other_points_dict[i] = int(order/gcd(other_points_dict[i] , order))
            return [order , other_points_dict]
        else:
            other_points_dict[tuple(vector_mod_Z(mat**-1*(order*point)))] = order

def order_of_new_equiv_point(point , dict_of_classes , new_equiv , order_of_equiv):
    zero_vec = zeros(len(new_equiv) , 1)
    order = 0
    other_points_dict = {}
    while True:
        order +=1
        for mult in range(order_of_equiv):
            res_vec = vector_mod_Z(order*point +mult*new_equiv)
            for other_vec in dict_of_classes[tuple(zero_vec)]:
                if res_vec == other_vec:
                    for i in other_points_dict:
                        other_points_dict[i] = int(order/gcd(other_points_dict[i] , order))
                    return [order , other_points_dict]

            for vec_class in dict_of_classes:
                for another_vec in vec_class:
                    if res_vec == another_vec:
                        other_points_dict[vec_class] = order

def classify_matrix(mat , list_of_points):
    [dict_of_points, max_order , max_ord_points , min_vec , dict_of_inv_points , min_vec_inv] = classify_points(mat , list_of_points)
    abelian_grps = [max_order]
    generators = [min_vec]
    total_ord = max_order
    if max_order == abs(det(mat)):
        gen_list = []
        for gen in generators:
            gen_list += (list(gen))
        B = Matrix(len(generators) , len(generators[0]), gen_list )

        Q = B*mat**-1*B.T
        return [abelian_grps , generators , Q]
    list_of_inv_points = []
    for key in dict_of_inv_points:
        list_of_inv_points.append(Matrix(list(key)))
    [dict_classes , classes , repr ] = quotient_space(list_of_inv_points , [min_vec_inv] , [max_order] )
    while len(dict_classes) >1:
        [dict_order_classes, max_order , max_ord_points , min_vec_inv] = classify_new_classes(dict_classes , repr , min_vec_inv , max_order )

        abelian_grps.append(max_order)
        total_ord *= max_order
        if dict_of_points[tuple(mat*min_vec_inv)] == max_order:
            generators.append(mat*min_vec_inv)
        else:
            for vec in dict_classes[tuple(min_vec_inv)]:
                if dict_of_points[tuple(mat*vec)] == max_order:
                    generators.append(mat*vec)
                    break
        if total_ord == det(mat):
            gen_list = []
            for gen in generators:
                gen_list += (list(gen))
            B = Matrix(len(generators) , len(generators[0]), gen_list )

            Q = B*mat**-1*B.T
            return [abelian_grps , generators , Q]
        [dict_classes , lst_new_classes , repr] = new_quot_space(dict_classes , repr , min_vec_inv , max_order)
    gen_list = []
    for gen in generators:
        gen_list += (list(gen))

    B = Matrix(len(generators) , len(generators[0]), gen_list )
    Q = B*mat**-1*B.T
    return [abelian_grps , generators , Q]

def compute_cokernel(mat):
    mat_sep = split_matrix(mat)
    mat_reduced = sub_matrix(mat_sep , mat.rank())
    set_of_points = set_of_points_within(mat_reduced)
    b = classify_matrix(mat_reduced , set_of_points)
    b.append(len(mat.row(0)) - mat.rank())
    b.append(det(mat_reduced))
    return b
