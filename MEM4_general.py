# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

from pqtrees.common_intervals.common_interval import CommonInterval
from pqtrees.pqtree_helpers.generate_s import IntervalHierarchy
from pqtrees.utilities.perm_helpers import tmap
from pqtrees.pqtree import PQTreeBuilder, PQTreeVisualizer
from pqtrees.common_intervals.preprocess_find import common_k_indexed_with_singletons
from pqtrees.pqtree_helpers.reduce_intervals import ReduceIntervals
from pqtrees.utilities.string_mutations import mutate_collection
from pqtrees.common_intervals.trivial import trivial_common_k_with_singletons
from copy import deepcopy

import glob

import math

import os


def get_pattern_from_csb_line(csb_line, csbs_names_dict):
    first_csbs_list = csb_line.split('\t')[4].split(',')
    csb_list = [csbs_names_dict[cog[:-1]] for cog in first_csbs_list]
    return tuple(csb_list)


def get_pattern_list_from_csb_line(csb_line, csbs_names_dict):
    first_csbs_list = csb_line.split('\t')[4].split(',')
    pattern = []
    for cog in first_csbs_list:
        pattern.append((str(csbs_names_dict[cog[:-1]]), str(cog[-1])))
    return pattern


def print_pqtree_from_csbs_file_path(csbs_file_path):
    csbs_file_lines = open(csbs_file_path).readlines()
    # #print(csbs_file_lines[1].split('\t')[4].split(','))
    first_csbs_list = csbs_file_lines[1].split('\t')[4].split(',')
    csbs_names_dict = {}
    for i in range(len(first_csbs_list)):
        csbs_names_dict[first_csbs_list[i][:-1]] = i

    perms = []
    for line in csbs_file_lines[1:]:
        perms.append(get_pattern_from_csb_line(line, csbs_names_dict))

    # #print(perms)

    strs = {"".join(str(x) for x in p) for p in perms}

    common_intervals_trivial = trivial_common_k_with_singletons(*perms)
    common_intervals = common_k_indexed_with_singletons(*perms)

    ir_intervals = ReduceIntervals.reduce(common_intervals)
    s = IntervalHierarchy.from_irreducible_intervals(ir_intervals)

    pqtree = PQTreeBuilder._from_s(s)
    # #print(pqtree.to_parens())


def get_pqtree_from_csbs_file_path_by_index(csbs_file_path, index):
    csbs_file_lines = open(csbs_file_path).readlines()
    first_csbs_list = csbs_file_lines[1].split('\t')[4].split(',')
    csbs_names_dict = {}
    for i in range(len(first_csbs_list)):
        csbs_names_dict[first_csbs_list[i][:-1]] = i

    perms = []
    for line in csbs_file_lines[1:]:
        perms.append(get_pattern_from_csb_line(line, csbs_names_dict))
    # #print(perms)

    # print(perms[0])

    perms, index_dict, index_dict_reverse = get_perms_by_index(perms, index)


    strs = {"".join(str(x) for x in p) for p in perms}

    common_intervals_trivial = trivial_common_k_with_singletons(*perms)
    common_intervals = common_k_indexed_with_singletons(*perms)

    ir_intervals = ReduceIntervals.reduce(common_intervals)
    s = IntervalHierarchy.from_irreducible_intervals(ir_intervals)

    pqtree = PQTreeBuilder._from_s(s)
    # print(pqtree.to_parens())
    #print(get_pqtree_by_dict(pqtree.to_parens(), index_dict_reverse))
    # PQTreeVisualizer.show(pqtree)
    return get_pqtree_by_dict(pqtree.to_parens(), index_dict_reverse)


def get_perms_by_index(perms, index):
    index_dict = {}
    index_dict_reverse = {}
    for i in range(len(perms[index])):
        index_dict[perms[index][i]] = i
        index_dict_reverse[i] = perms[index][i]

    new_perms = [tuple([index_dict[cog] for cog in perms[index]])]
    for i in range(len(perms)):
        if i!=index:
            new_perms.append(tuple([index_dict[cog] for cog in perms[i]]))

    return new_perms, index_dict, index_dict_reverse


def get_pqtree_by_dict(pqtree, index_dict):
    to_return = ""
    for c in pqtree:
        if c!='[' and c!=']' and c!='(' and c!=')' and c!=' ':
            to_return = to_return + str(index_dict[int(c)])
        else:
            to_return = to_return + c
    return to_return


class PQNode:
    def __init__(self, type, children,span, direction, flipped, l, r, level):
        self.type = type
        self.children = children
        self.span = span
        self.direction = direction
        self.flipped = flipped
        self.l = l
        self.r = r
        self.level = level

    def __str__(self):
        to_str = ""
        for child in self.children:
            to_str = to_str +" "+ str(child)

        if self.type == "P":
            to_str = "(" + to_str[1:] + ")"
        else:
            to_str = "[" + to_str[1:] + "]"
        return to_str


class Leaf:
    def __init__(self, cog, direction, flipped, l, r):
        self.cog = cog
        self.direction = direction
        self.flipped = flipped
        self.l = l
        self.r = r
        self.span = 1

    def __str__(self):
        if self.direction == 1:
            return str(self.cog) + "+"
        return str(self.cog) + "-"

def get_direction(children_list):
    sum_pos = sum([child.direction==1 for child in children_list])
    if sum_pos >= len(children_list):
        return 1
    return -1


def get_level_by_children(children_list):
    highest_level = 1
    for child in children_list:
        if type(child) is PQNode:
            if child.level + 1 > highest_level:
                highest_level = child.level + 1

    return highest_level


def get_pqtree_by_parens(parens, directions):
    # print("PPPPPPPPPPPPPPPPPPPPPPPPPPPP",parens)
    p_stack = [")"]
    q_stack = ["]"]
    children_list = []
    # parens = parens.replace(" ", "")
    # #print(parens)
    span = 0
    for i in range(1, len(parens)-1):
        if parens[i]=="(" and q_stack[-1]=="]":
            p_stack.append(i)
        if parens[i] == "[" and p_stack[-1]==")":
            q_stack.append(i)
        if parens[i] == ")" and q_stack[-1]=="]":
            child_parens = parens[p_stack.pop():i+1]
            if p_stack[-1]==")":
                child_node = get_pqtree_by_parens(child_parens, directions)
                span = span + child_node.span
                children_list.append(child_node)
        if parens[i] == "]" and p_stack[-1]==")":
            child_parens = parens[q_stack.pop():i+1]
            if q_stack[-1]=="]":
                child_node = get_pqtree_by_parens(child_parens, directions)
                span = span + child_node.span
                children_list.append(child_node)
        if parens[i] != "(" and parens[i] != ")" and parens[i] != "[" and parens[i] != "]" and p_stack[-1]==")" and q_stack[-1]=="]":
            span = span + 1
            children_list.append(Leaf(parens[i], directions[parens[i]], False, 0, 0))

    direction = get_direction(children_list)
    level = get_level_by_children(children_list)
    if parens[0]=="(":
        return PQNode("P", children_list, span, direction, False, 0, 0, level)
    else:
        # #   newwwwww - if there is a Qnode that has 2 children and its in high level, it should be Pnode
        # if len(children_list) == 2 and span != 2:
        #     return PQNode("P", children_list, span, direction, False, 0, 0, level)
        # #   newwwwww
        return PQNode("Q", children_list, span, direction, False, 0, 0, level)


def get_pattern_by_file_and_index(csbs_file_path, index):
    csbs_file_lines = open(csbs_file_path).readlines()
    # #print(csbs_file_lines[1].split('\t')[4].split(','))
    first_csbs_list = csbs_file_lines[1].split('\t')[4].split(',')
    csbs_names_dict = {}
    for i in range(len(first_csbs_list)):
        csbs_names_dict[first_csbs_list[i][:-1]] = i

    line = csbs_file_lines[index+1]
    return get_pattern_list_from_csb_line(line, csbs_names_dict)


def get_pattern_directions_by_file_and_index(csbs_file_path, index):
    csbs_file_lines = open(csbs_file_path).readlines()
    # #print(csbs_file_lines[1].split('\t')[4].split(','))
    first_csbs_list = csbs_file_lines[1].split('\t')[4].split(',')
    csbs_names_dict = {}
    for i in range(len(first_csbs_list)):
        csbs_names_dict[first_csbs_list[i][:-1]] = i

    line = csbs_file_lines[index+1]

    dict = {}
    csb_splited = line.split('\t')[4].split(',')
    for cog in csb_splited:
        if cog[-1] == "+":
            dict[str(csbs_names_dict[cog[:-1]])] = 1
        else:
            dict[str(csbs_names_dict[cog[:-1]])] = -1

    return dict


def get_pqtree_from_file_and_index(path, index):
    pq_tree_parens = get_pqtree_from_csbs_file_path_by_index(path, index)
    pq_tree_parens = pq_tree_parens.replace(" ", "")
    csb_directions = get_pattern_directions_by_file_and_index(path, index)
    pq_tree = get_pqtree_by_parens(pq_tree_parens, csb_directions)
    return pq_tree


# function to generate all the sub lists
def sub_lists(l):
    lists = [[]]
    for i in range(len(l) + 1):
        for j in range(i):
            lists.append(l[j: i])
    return lists


def is_vertex_cover(vertex_set, edges):
    for edge in edges:
        if edge[0] not in vertex_set and edge[1] not in vertex_set:
            return False
    return True


def get_vertex_cover_score(vertex_cover, weighted_vertex_dict):
    sum = 0
    for vertex in vertex_cover:
        sum = sum + weighted_vertex_dict[vertex]
    return sum


def get_minimum_weighted_vertex_cover(vertexes, edges, weighted_vertex_dict):
    min_value = float('inf')
    vertex_cover = 0
    all_sub_lists = sub_lists(vertexes)
    for sub_list in all_sub_lists:
        if is_vertex_cover(set(sub_list), edges):
            score = get_vertex_cover_score(sub_list, weighted_vertex_dict)
            if score < min_value:
                vertex_cover = sub_list
                min_value = score

    return min_value, vertex_cover





















from itertools import chain, combinations

def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    ps = list(chain.from_iterable(frozenset(combinations(s, r)) for r in range(len(s)+1)))
    ps = [frozenset(e) for e in ps]
    return ps








def get_bp_from_children_orders(co_tree, co_string):
  # #print(co_tree, co_string)
  bp = 0
  for i in range(1,len(co_tree)):
    tup1 = 0
    if co_tree[i-1][1] == 1 and co_tree[i][1] == 1:
        tup1 = (co_tree[i - 1][0], co_tree[i][0])
    if co_tree[i - 1][1] == -1 and co_tree[i][1] == -1:
        tup1 = (co_tree[i][0], co_tree[i - 1][0])
    if tup1 == 0:   # one is flipped and one is not: its a breakpoint
        # #print("hereeeeeeeeeeeeee")
        bp += 1
        continue
    found = False
    for j in range(1,len(co_string)):
      tup2 = (co_string[j-1], co_string[j])
      if tup1 == tup2:
        found = True

    if(not found):

      bp += 1

  # #print(bp)
  return bp


def get_jumping_penalty(co_tree, co_string, children_span_list):
    vertexes = [i for i in range(len(co_string))]
    edges = []
    for i in range(0, len(co_tree)-1):
        for t in range(i+1, len(co_tree)):
            tup1 = 0
            if co_tree[i][1] == 1 and co_tree[t][1] == 1:
                tup1 = (co_tree[i][0], co_tree[t][0])
            if co_tree[i][1] == -1 and co_tree[t][1] == -1:
                tup1 = (co_tree[t][0], co_tree[i][0])
            if tup1 == 0:  # one is flipped and one is not: its a breakpoint
                continue

            is_jumped = False
            for cog in co_string:
                if cog == tup1[1]:
                    is_jumped = True
                    break
                if cog == tup1[0]:
                    break

            if is_jumped:
                edges.append((int(tup1[0]), int(tup1[1])))

    weighted_vertex_list = [(span-1)/2 for span in children_span_list]
    minimum_VC, a = get_minimum_weighted_vertex_cover(vertexes, edges, weighted_vertex_list)

    return minimum_VC


def check_if_all_children_flipped(node):
    for child in node.children:    # checking if all children were flipped
        if not child.flipped:
            return False
    return True


def check_if_all_children_changed_direction(node):
    for child in node.children:    # checking if all children were changed directions
        if type(child) is Leaf:
            if not child.flipped:
                return False
        else:
            if child.direction == 1:
                return False
    return True


def get_children_flipped_penalties(C):
    penalties_sum = 0
    for child in C:
        if type(child) is Leaf:
            # penalties_sum += qnode_flip_penalty*1
            penalties_sum += qnode_flip_penalty
        else:
            if child.type == "Q":
                # penalties_sum += qnode_flip_penalty + math.log(child.span)
                penalties_sum += qnode_flip_penalty
    return penalties_sum


def check_if_node_children_flip_penalties(pnode, i, e):

    index = i
    flipped = True
    for j in range(len(pnode.children) - 1, -1, -1): # checking children from right to left
        child = pnode.children[j]
        if child.l == index:
            index = child.r + 1
        else:
            flipped = False
    if flipped == False or index > e+1:
        flipped = False

    if flipped == False:
        return False



    all_children_flipped = check_if_all_children_flipped(pnode)
    all_children_changed_direction = check_if_all_children_changed_direction(pnode)

    if flipped and all_children_changed_direction and all_children_flipped:
        return True

    return False


def is_VC_cover_y(C, C_tag, N, y, node):
    if y in C_tag:
        return True
    positions = {node.children[i]: i for i in range(len(node.children))}
    changed_direction = {child: ((child.direction==-1 and child not in N) or (child.direction==1 and child in N)) for child in node.children}
    for x in C:
        if x not in C_tag:
            if positions[x] < positions[y] and changed_direction[x] and changed_direction[y]:
                return False
            if positions[x] > positions[y] and not changed_direction[x] and not changed_direction[y]:
                return False

    return True


def is_valid_entry(C, C_tag, N, k_T, k_S, y, pqnode, i, e, pos, csb, end_point):
    # print("csbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb")
    # print(csb)
    if end_point > len(csb) - 1:
        return False
    if end_point > e:
        return False
    if length_L(C,k_T,0) < 0:   # K_T > span(C)
        return False
    length_of_derivation = length_L(C,k_T,k_S)
    if k_S > length_of_derivation or length_of_derivation > e - i + 1:
        return False
    if not is_VC_cover_y(C, C_tag, N, y, pqnode):
        return False
    if len(C) == 1 and (((pos > span({y})/2) and (y in N)) or ((pos < span({y})/2) and (y not in N))):
        return False
    # if ((pos > span({y})/2) and (y in N)) or ((pos < span({y})/2) and (y not in N)):
    #     return False
    return True


def is_valid_entry_Q(C, N, k_T, k_S, pqnode, i, e, pos, xj, csb, end_point):
    if end_point > len(csb) - 1:
        return False
    if length_L(C,k_T,0) < 0:   # K_T > span(C)
        return False
    length_of_derivation = length_L(C,k_T,k_S)
    if k_S > length_of_derivation or length_of_derivation > e - i + 1:
        return False
    if len(C) == 1 and (((pos > span({xj})/2) and (xj in N)) or ((pos < span({xj})/2) and (xj not in N))):
        return False
    # if ((pos > span({xj})/2) and (xj in N)) or ((pos < span({xj})/2) and (xj not in N)):
    #     return False
    return True


def check_is_valid_entry(C, C_tag, N, k_T, k_S, y, pqnode, i, e, pos, csb):
    if length_L(C,k_T,0) < 0:   # K_T > span(C)
        print("11111111111111111111111111111111111111111111111111111111111111111")
    length_of_derivation = length_L(C,k_T,k_S)
    if k_S > length_of_derivation or length_of_derivation > e - i + 1:
        print("22222222222222222222222222222222222222222222222222222222222222222")
    if not is_VC_cover_y(C, C_tag, N, y, pqnode):
        print("33333333333333333333333333333333333333333333333333333333333333333")
    if len(C) == 1 and (((pos > span({y})/2) and (y in N)) or ((pos < span({y})/2) and (y not in N))):
        print("44444444444444444444444444444444444444444444444444444444444444444", pos)



def get_y_derivations_dict(y, y_sign, end_point, k_T, k_S, sos_pos, A):
    derivations_dict = {}
    for r_T in range(k_T+1):
        max_pos = span({y}) - r_T
        for r_S in range(k_S+1):
            for pos in range(max_pos + 1):
                y_start_index = end_point - length_L({y}, r_T, r_S) + 1
                if length_L({y}, r_T, r_S) > 0:     # the derivation is not empty
                    if type(y) is Leaf:
                        if (y_sign < 0 and pos == 0) or (y_sign > 0 and pos == 1):
                            # print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX", end_point, length_L({y}, r_T, r_S), y_start_index)
                            derivations_dict[(r_T, r_S, pos)] = A[y][y_start_index][r_T][r_S][pos]
                    else:
                        if (y_sign < 0 and pos <= max_pos/2) or (y_sign > 0 and pos >= max_pos/2):
                            derivations_dict[(r_T, r_S, pos)] = A[y][y_start_index][r_T][r_S][pos]

    return derivations_dict


def dist_delta(N, y, z, node):  # for case 2
    positions = {node.children[i]: i for i in range(len(node.children))}
    changed_direction = {child: ((child.direction==-1 and child not in N) or (child.direction==1 and child in N)) for child in node.children}
    if positions[z]+1 == positions[y] and not changed_direction[z] and not changed_direction[y]:
        return 0
    if positions[z]-1 == positions[y] and changed_direction[z] and changed_direction[y]:
        return 0

    # print("distdeltaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa " + str(y.cog) + str(y.direction) + " " + str(z.cog) + str(z.direction))
    # return 1 / node.level
    return 1


def dist_delta2(N, y, z, node):  # for case 3
    positions = {node.children[i]: i for i in range(len(node.children))}
    changed_direction = {child: ((child.direction==-1 and child not in N) or (child.direction==1 and child in N)) for child in node.children}
    if positions[z] < positions[y] and not changed_direction[z] and not changed_direction[y]:
        return 0
    if positions[z] > positions[y] and changed_direction[z] and changed_direction[y]:
        return 0

    return 1


def dist_delta_Q(N, xj, xj_minus1, is_l):
    if is_l:
        changed_direction = {
            child: ((child.direction == -1 and child not in N) or (child.direction == 1 and child in N)) for child in
            [xj, xj_minus1]}
    else:
        changed_direction = {child: ((child.direction == 1 and child not in N) or (child.direction == -1 and child in N)) for child in [xj, xj_minus1]}

    if not changed_direction[xj] and not changed_direction[xj_minus1]:
        return 0

    return 1

def dist_delta_Q2(N, xj, xr, is_l):
    if is_l:
        changed_direction = {
            child: ((child.direction == -1 and child not in N) or (child.direction == 1 and child in N)) for child in
            [xj, xr]}
    else:
        changed_direction = {child: ((child.direction == 1 and child not in N) or (child.direction == -1 and child in N)) for child in [xj, xr]}

    if not changed_direction[xj] and not changed_direction[xr]:
        return 0

    return 1

def jump_violation_delta(C_tag,y):
    if y in C_tag:
        return jump_penalty * ((y.span-1)/2)
    return 0


def get_Cyz_dict(C,y,node):
    positions = {node.children[i]: i for i in range(len(node.children))}
    Cyz_dict = {}
    for z in C:
        if z != y:
            min_pos, max_pos = min([positions[z],positions[y]]), max([positions[z],positions[y]])
            Cyz = frozenset(node.children[min_pos+1:max_pos])
            if Cyz.issubset(C):
                Cyz_dict[z] = Cyz
    return Cyz_dict


def span(U):    # U needs to be iterable
    return sum([node.span for node in U])


def get_P(P, C, C_tag, N, k_T, k_S, y, pos):
    if k_T < 0 or k_S < 0 or pos < 0 or pos > span(C) - k_T:
        return float('inf')
    return P[C][C_tag][N][k_T][k_S][y][pos]


def get_Q(Q, span_C, j, k_T, k_S, N, pos):
    if k_T < 0 or k_S < 0 or pos < 0 or pos > span_C - k_T:
        return float('inf')
    if j == 0:
        return 0
    return Q[j][k_T][k_S][N][pos]



def check_dist_delta(N, y, pnode, C_minus_y):
    for z in C_minus_y:
        if dist_delta(N, y, z, pnode)==0:
            print("00000000000000000000000000000000000000000000000000000000000   " + str(y.cog) + "   " + str(z.cog))
            print(z.cog)
        if dist_delta(N, y, z, pnode)==1:
            print("111111111111111111111111111111111111111111111111111111111111   " + str(y.cog) + "   " + str(z.cog))
            print(z.cog)


def p_mapping(pnode, i, e, csb, d_T, d_S, A):
    P = {}
    # print("NOTVALIDDDDDDDDDDDDDDDDDDDDDDDDDD")
    for C in powerset(pnode.children):
        P[C] = {}
        for C_tag in powerset(C):
            P[C][C_tag] = {}
            for N in powerset(C):
                P[C][C_tag][N] = {}
                for k_T in range(d_T+1):
                    P[C][C_tag][N][k_T] = {}
                    for k_S in range(d_S+1):
                        P[C][C_tag][N][k_T][k_S] = {}
                        for y in C:
                            P[C][C_tag][N][k_T][k_S][y] = {}
                            for sos_pos in range(span(C) + 1 - k_T):
                                C_minus_y, C_tag_minus_y, N_minus_y = C.difference({y}), C_tag.difference({y}), N.difference({y})
                                y_sign = (y not in N)*2 - 1
                                end_point = end_point_E(C, i, k_T, k_S)
                                # check_dist_delta(N, y, pnode, C_minus_y)
                                if is_valid_entry(C, C_tag, N, k_T, k_S, y, pnode, i, e, sos_pos, csb, end_point):
                                    # if pnode.children[0] in N and pnode.children[1] in N and pnode.children[2] not in N and len(C)==3:
                                    #     print("innnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn  " + str(dist_delta(N, pnode.children[0], pnode.children[1], pnode)))
                                    # if len(C_tag)<1:
                                    #     print("OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO")
                                    if len(C) == 1:
                                        y_start_index = end_point - length_L({y}, k_T, k_S) + 1
                                        P[C][C_tag][N][k_T][k_S][y][sos_pos] = A[y][y_start_index][k_T][k_S][sos_pos] + jump_violation_delta(C_tag,y)
                                    else:
                                        y_derivations_dict = get_y_derivations_dict(y, y_sign, end_point, k_T, k_S, sos_pos, A)

                                        # case1 = P[C][C_tag][N][k_T][k_S-1][y][sos_pos] + delete_S_penalty
                                        case1 = get_P(P, C, C_tag, N, k_T, k_S-1, y, sos_pos) + delete_S_penalty
                                        # case2 = min([P[C_minus_y][C_tag_minus_y][N_minus_y][k_T-deriv[0]][k_S-deriv[1]][z][sos_pos-deriv[2]] + score + dist_delta(N, y, z, pnode) + jump_violation_delta(C_tag,y) for (deriv, score) in y_derivations_dict.items() for z in C_minus_y])
                                        case2 = float('inf')
                                        if len(y_derivations_dict.items()) > 0:
                                            case2 = min([get_P(P, C_minus_y, C_tag_minus_y, N_minus_y, k_T-deriv[0], k_S-deriv[1], z, sos_pos-deriv[2]) + score + dist_delta(N, y, z, pnode) + jump_violation_delta(C_tag,y) for (deriv, score) in y_derivations_dict.items() for z in C_minus_y])
                                        # case3 = min([P[C_minus_y.difference(Cyz)][C_tag_minus_y.difference(Cyz)][N_minus_y.difference(Cyz)][k_T-deriv[0]-span(Cyz)][k_S-deriv[1]][z][sos_pos-deriv[2]] + score + dist_delta2(N, y, z, pnode) + jump_violation_delta(C_tag,y) for (deriv, score) in y_derivations_dict.items() for (z,Cyz) in get_Cyz_dict(C,y,pnode).items()])
                                        case3 = float('inf')
                                        # if len(get_Cyz_dict(C,y,pnode).items()) > 0:
                                        #     case3 = min([get_P(P, C_minus_y.difference(Cyz), C_tag_minus_y.difference(Cyz), N_minus_y.difference(Cyz), k_T-deriv[0]-span(Cyz), k_S-deriv[1], z, sos_pos-deriv[2]) + score + dist_delta2(N, y, z, pnode) + jump_violation_delta(C_tag,y) for (deriv, score) in y_derivations_dict.items() for (z,Cyz) in get_Cyz_dict(C,y,pnode).items()])
                                        P[C][C_tag][N][k_T][k_S][y][sos_pos] = min([case1, case2, case3])

                                        # if i == 6 and type(y) is Leaf and y.cog == "3":
                                        #     if len(C) == 2 and pnode.children[1] in C and pnode.children[1].cog == "4" and pnode.children[2].cog == "5" and pnode.children[1] in N and y not in N:
                                        #         print("INNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN", dist_delta(N, y, pnode.children[1], pnode), "C_tag: ", len(C_tag))

                                        # if y == pnode.children[0] and pnode.children[0] in N and pnode.children[1] in N and pnode.children[2] not in N and len(C) == 3:
                                        #     print("innnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn  " + str(dist_delta(N, pnode.children[0], pnode.children[1], pnode)) + "   " + str(P[C][C_tag][N][k_T][k_S][y][sos_pos])+ "   " + str(P[C.difference({pnode.children[0]})][C_tag.difference({pnode.children[0]})][N.difference({pnode.children[0]})][k_T][k_S][pnode.children[1]][1]))

                                        # if y == pnode.children[1] and pnode.children[2] in C and pnode.children[0] not in C and y in N and pnode.children[2] not in N and sos_pos==1:
                                        #     print("innnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn  " + str(P[C][C_tag][N][k_T][k_S][y][sos_pos]) + "  " + str(y_derivations_dict[(0, 0, 0)]))

                                        # if y == pnode.children[1] and pnode.children[2] in C and pnode.children[0] not in C and y in N and pnode.children[2] not in N:
                                        #     print("innnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn  " + str(P[C][C_tag][N][k_T][k_S][y][sos_pos]))

                                        # if i == 0 and type(y) is PQNode and span({pnode})==9 and y == pnode.children[0] and pnode.children[2] in C and len(C) == 2 and len(N) == 0:
                                        #     print("INNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN  ", len(C_tag), jump_violation_delta(C_tag,y))


                                else:
                                    # if i == 0 and type(y) is PQNode and span({pnode}) == 9 and y == pnode.children[2] and len(C) == 1 and len(N) == 0:
                                    #     print("INNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN  ", len(C_tag), jump_violation_delta(C_tag, y))
                                    #     check_is_valid_entry(C, C_tag, N, k_T, k_S, y, pnode, i, e, sos_pos, csb)
                                    P[C][C_tag][N][k_T][k_S][y][sos_pos] = float('inf')

    for k_T in range(d_T+1):
        for k_S in range(d_S+1):
            for pos in range(span({pnode}) + 1 - k_T):
                #TODO: min_val...... like in q-mapping
                min_val = float('inf')
                for C in powerset(pnode.children):
                    if span({pnode}) - span(C) <= k_T:
                        min_val_C = min([get_P(P, C, C_tag, N, k_T, k_S, y, pos) + delete_T_penalty*span(set(pnode.children).difference(C)) for C_tag in powerset(C) for N in powerset(C) for y in C])
                        min_val = min([min_val, min_val_C])
                        # print("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA  " + str(A[pnode][i][k_T][k_S][pos]) + "  " + str(pos) + "  " + str(len(C)))
                        # checkentry(P, C, pos, pnode)

                A[pnode][i][k_T][k_S][pos] = min_val
    return A


def checkentry(P, C, pos, pnode):
    i=0
    j=0
    k=0
    for C_tag in powerset(C):
        i += 1
        j=0
        k=0
        for N in powerset(C):
            j += 1
            k = 0
            for y in C:
                k += 1
                if get_P(P, C, C_tag, N, 0, 0, y, pos) == 0 and pnode.children[0]==y:
                    print("CHECKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK", i, j, k, len(N), pos)


def q_mapping(qnode, i, e, csb, d_T, d_S, A):
    Ql = {}
    Qr = {}
    for j in range(1, len(qnode.children)+1):
        Ql[j], Qr[j] = {}, {}
        Cl, Cr = qnode.children[:j], qnode.children[-j:]
        Cr.reverse()

        Ql = q_iteration(qnode, i, e, csb, d_T, d_S, Ql, Cl, j, A, True)

        Qr = q_iteration(qnode, i, e, csb, d_T, d_S, Qr, Cr, j, A, False)

    for k_T in range(d_T+1):
        for k_S in range(d_S+1):
            for pos in range(span({qnode}) + 1 - k_T):
                min_val = float('inf')
                for j in range(1, len(qnode.children)+1):
                    Cl, Cr = qnode.children[:j], qnode.children[-j:]
                    Cr.reverse()
                    span_Cl, span_Cr = span(Cl), span(Cr)
                    span_d_Cl, span_d_Cr = qnode.span-span_Cl, qnode.span-span_Cr
                    min_l, min_r = float('inf'), float('inf')

                    if span({qnode}) - span_Cl <= k_T:
                        min_l = min([get_Q(Ql, span_Cl, j, k_T-span_d_Cl, k_S, N, pos) + delete_T_penalty*span_d_Cl for N in powerset(Cl)])

                    if span({qnode}) - span_Cr <= k_T:
                        min_r = min([get_Q(Qr, span_Cr, j, k_T-span_d_Cr, k_S, N, pos) + delete_T_penalty*span_d_Cr + flip_check(Cr, N) for N in powerset(Cr)])

                    min_val = min([min_val, min_l, min_r])
                    # if span(Cl) == 3:
                    #     print("QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ", min_l, min_r, "   ", get_Q(Qr, span_Cr, j, k_T-span_d_Cr, k_S, frozenset(Cr), pos), pos)

                    # if j==len(qnode.children) and e==1:
                    #     print("MINVALLLLLLLLLLLLLLLLLLLLLLLLLLLLLL", min_l, min_r)

                A[qnode][i][k_T][k_S][pos] = min_val

    return A


def q_iteration(qnode, i, e, csb, d_T, d_S, Q, C, j, A, is_l):
    # span_C = span(C)
    span_C = {r: span([c for c in C[:r]]) for r in range(j+1)}
    span_r_to_j = {r: span([c for c in C[r+1:j-1]]) for r in range(j)}
    for k_T in range(d_T + 1):
        Q[j][k_T] = {}
        for k_S in range(d_S + 1):
            Q[j][k_T][k_S] = {}
            for N in powerset(C):
                Q[j][k_T][k_S][N] = {}
                for pos in range(span(C) + 1 - k_T):
                    xj = C[j-1]
                    N_minus_xj = N.difference({xj})
                    N_r = {r: N.intersection(frozenset(C[:r])) for r in range(j)}
                    # N_r = {r: N.difference(frozenset(C[:j])) for r in range(j)}
                    xj_sign = (xj not in N) * 2 - 1
                    xj_span = xj.span
                    end_point = end_point_E(C, i, k_T, k_S)
                    if is_valid_entry_Q(C, N, k_T, k_S, qnode, i, e, pos, xj, csb, end_point):
                        if j == 1:
                            xj_start_index = end_point - length_L({xj}, k_T, k_S) + 1
                            Q[j][k_T][k_S][N][pos] = A[xj][xj_start_index][k_T][k_S][pos]
                            # if not is_l:
                            #     print("XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX",xj, A[xj][xj_start_index][k_T][k_S][pos])
                        else:
                            xj_derivations_dict = get_y_derivations_dict(xj, xj_sign, end_point, k_T, k_S, pos, A)
                            xj_minus1 = C[j - 2]

                            case1 = get_Q(Q, span_C[j], j, k_T, k_S-1, N, pos) + delete_S_penalty

                            # case2 = get_Q(Q, span_C[j]-xj_span, j-1, k_T-xj_span, k_S, N_minus_xj, pos) + xj_span*delete_T_penalty
                            #
                            # case3 = float('inf')
                            # if len(xj_derivations_dict.items()) > 0:
                            #     case3 = min([get_Q(Q, span_C[j]-xj_span, j-1, k_T - deriv[0], k_S - deriv[1], N_minus_xj, pos - deriv[2]) + score + bp_qnode_penalty*dist_delta_Q(N, xj, xj_minus1, is_l) for (deriv, score) in xj_derivations_dict.items()])


                            case2 = float('inf')
                            if len(xj_derivations_dict.items()) > 0:
                                case2 = min([get_Q(Q, span_C[r], r, k_T - deriv[0] - span_r_to_j[r], k_S - deriv[1], N_r[r], pos - deriv[2]) + score + bp_qnode_penalty*dist_delta_Q2(N, xj, C[r-1], is_l) + delete_T_penalty*span_r_to_j[r] for r in [j-1] for (deriv, score) in xj_derivations_dict.items()]) # no deletions
                                # case2 = min([get_Q(Q, span_C[r], r, k_T - deriv[0] - span_r_to_j[r], k_S - deriv[1], N_r[r], pos - deriv[2]) + score + bp_qnode_penalty*dist_delta_Q2(N, xj, C[r-1], is_l) + delete_T_penalty*span_r_to_j[r] for r in range(1, j) for (deriv, score) in xj_derivations_dict.items()])

                            case3 = float('inf')

                            Q[j][k_T][k_S][N][pos] = min([case1, case2, case3])


                            # if j == 3 and xj in N and xj_minus1 in N and len(N)==3 and not is_l:
                            #     xj_start_index = end_point - length_L({xj}, k_T, k_S) + 1
                            #     print("YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY", xj, A[xj][xj_start_index][k_T][k_S][0], A[xj][xj_start_index][k_T][k_S][1], "  ", Q[j][k_T][k_S][N][pos])

                    else:
                        Q[j][k_T][k_S][N][pos] = float('inf')
    return Q


def flip_check(Cr, N):
    changed_direction = {child: ((child.direction == 1 and child not in N) or (child.direction == -1 and child in N))
                         for child in Cr}
    for child in Cr:
        if changed_direction[child]:
            return 0

    children_flipped_penalties = get_children_flipped_penalties(Cr)
    # flip_penalty = qnode_flip_penalty + math.log(span(Cr))
    flip_penalty = qnode_flip_penalty
    # print("HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH", children_flipped_penalties, flip_penalty, span(Cr))
    # print("FLIPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP", span(Cr), flip_penalty - children_flipped_penalties)
    return flip_penalty - children_flipped_penalties
    # return 0


def length_L(nodes, k_T, k_S):
    sum_of_spans = sum([node.span for node in nodes])
    return sum_of_spans - k_T + k_S


def end_point_E(nodes, i, k_T, k_S):
    return i - 1 + length_L(nodes, k_T, k_S)



def iteration_leaf(leaf, csb, d_T, d_S, A):
    n = len(csb)
    for i in range(n):
        for k_S in range(d_S+1):
            for pos in [0,1]:
                if d_T > 0:
                    A[leaf][i][1][k_S][pos] = delete_T_penalty + k_S*delete_S_penalty    # the leaf and all k_S characters are deleted
                if (pos == 1 and csb[i][1] == "+") or (pos == 0 and csb[i][1] == "-"):
                    for j in range(i, min(i+k_S+1, len(csb))):
                        if csb[j][0] == leaf.cog: # TODO
                            A[leaf][i][0][k_S][pos] = k_S*delete_S_penalty   # all k_S characters are deleted
                            if (leaf.direction == 1 and pos == 0) or (leaf.direction == -1 and pos == 1):
                                A[leaf][i][0][k_S][pos] = A[leaf][i][0][k_S][pos] + qnode_flip_penalty

    # print("min val: " + str(min([A[leaf][i][0][k_S][pos] for i in range(n) for k_S in range(d_S+1) for pos in [0,1]])))
    return A


# def iteration_leaf(leaf, csb, d_T, d_S, A):
#     n = len(csb)
#     for i in range(n):
#         for k_S in range(d_S+1):
#             for pos in [0,1]:
#                 if d_T > 0:
#                     A[leaf][i][1][k_S][pos] = delete_T_penalty + k_S*delete_S_penalty    # the leaf and all k_S characters are deleted
#                 if (pos == 1 and csb[i][1] == "+") or (pos == 0 and csb[i][1] == "-"):
#                     for j in range(i, min(i+k_S+1, len(csb))):
#                         if csb[j][0] == leaf.cog: # TODO
#                             A[leaf][i][0][k_S][pos] = k_S*delete_S_penalty   # all k_S characters are deleted
#
#     # print("min val: " + str(min([A[leaf][i][0][k_S][pos] for i in range(n) for k_S in range(d_S+1) for pos in [0,1]])))
#     return A





def iteration_internal_node(pqtree_node, csb, d_T, d_S, A):
    for i in range(len(csb)):
        e = end_point_E([pqtree_node], i, 0, d_S)
        # if e > len(csb)-1:
        #     print("NOTVALIDDDDDDDDDDDDDDDDDDDDDDDDDD")
        #     continue
        if pqtree_node.type == "P":
            A = p_mapping(pqtree_node, i, e, csb, d_T, d_S, A)
        else:
            A = q_mapping(pqtree_node, i, e, csb, d_T, d_S, A)
    return A


def create_entries_A(pqtree_node, csb, d_T, d_S, A):
    A[pqtree_node] = {}
    for i in range(len(csb)):
        A[pqtree_node][i] = {}
        for k_T in range(d_T+1):
            A[pqtree_node][i][k_T] = {}
            for k_S in range(d_S + 1):
                A[pqtree_node][i][k_T][k_S] = {}
                for sos_pos in range(span({pqtree_node}) + 1 - k_T):
                    A[pqtree_node][i][k_T][k_S][sos_pos] = float('inf')
    return A


def main_algorithm(pqtree_node, csb, d_T, d_S, A):
    A = create_entries_A(pqtree_node, csb, d_T, d_S, A)
    if type(pqtree_node) is Leaf:
        A = iteration_leaf(pqtree_node, csb, d_T, d_S, A)
        # print("min val: " + str(min([A[pqtree_node][i][0][k_S][pos] for i in range(len(csb)) for k_S in range(d_S+1) for pos in [0,1]])))
    else:
        for child in pqtree_node.children:
            A = main_algorithm(child, csb, d_T, d_S, A)
        A = iteration_internal_node(pqtree_node, csb, d_T, d_S, A)
        # print("---------------------------------------------min val: " + str(min([A[pqtree_node][i][0][k_S][pos] for i in range(len(csb)) for k_T in range(d_T+1) for k_S in range(d_S + 1) for pos in range(span({pqtree_node}) + 1 - k_T)])))

    return A


def remove_and_get_line(fileName, lineToSkip):
    """ Removes a given line from a file """
    with open(fileName, 'r') as read_file:
        lines = read_file.readlines()

    line_skipped = 0
    currentLine = 0
    with open(fileName, 'w') as write_file:
        for line in lines:
            if currentLine == lineToSkip+1:
                line_skipped = line
            else:
                write_file.write(line)

            currentLine += 1

    return line_skipped


def write_line_at_index(fileName, index, lineToWrite):
    """ Removes a given line from a file """
    with open(fileName, 'r') as read_file:
        lines = read_file.readlines()


    currentLine = 0
    with open(fileName, 'w') as write_file:
        for line in lines:
            if currentLine == index:
                write_file.write(line)
                write_file.write(lineToWrite)
            else:
                write_file.write(line)

            currentLine += 1


def count_num_of_internal_nodes(pq_tree_string):
    sum = 0
    for c in pq_tree_string:
        if c == '[':
            sum += 1
    return sum


def get_num_of_csbs_in_file(path):
    num = 0
    with open(path, 'r') as read_file:
        lines = read_file.readlines()
        num = len(lines) - 1
    return num


def find_outliers(path):    # returns the outlier's index or -1 if not exist
    num_of_csbs = get_num_of_csbs_in_file(path)
    all_pq_tree_string = get_pqtree_from_csbs_file_path_by_index(path, 0)
    all_internal_nodes = count_num_of_internal_nodes(all_pq_tree_string)
    # #print("------------ all ", all_pq_tree_string, all_internal_nodes)

    outlier = -1
    num_outliers = 0
    for i in range(num_of_csbs):
        line_removed = remove_and_get_line(path, i)     # removing line at index i

        pq_tree_string = get_pqtree_from_csbs_file_path_by_index(path, 0)
        internal_nodes = count_num_of_internal_nodes(pq_tree_string)

        # #print("------------", i, pq_tree_string, internal_nodes)
        if internal_nodes > all_internal_nodes:
            outlier = i
            num_outliers += 1

        write_line_at_index(path, i, line_removed)     # writing the line at index i

    #print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! num of outliers:", num_outliers, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")

    # if num_outliers > 1:
    #     outlier = -1

    return outlier

def calculate_s_score_node(node):
    pnodes_childs = []
    qnodes_num = 0

    if node.type == "Q":
        qnodes_num = 1

    if node.type == "P":
        pnodes_childs = [len(node.children)]

    for child in node.children:
        if type(child) is Leaf:
            continue
        pnodes_childs_child, qnodes_num_cild = calculate_s_score_node(child)
        pnodes_childs.extend(pnodes_childs_child)
        qnodes_num += qnodes_num_cild


    return pnodes_childs, qnodes_num


def calculate_s_score(pq_tree):
    pnodes_childs, qnodes_num = calculate_s_score_node(pq_tree)
    pnodes_factorial_mult = 1
    for pnodes_num in pnodes_childs:
        pnodes_factorial_mult = pnodes_factorial_mult*math.factorial(pnodes_num)

    score_down = math.pow(2, qnodes_num)*pnodes_factorial_mult
    score_up = math.factorial(pq_tree.span)
    score = score_up/score_down
    return score

def find_outliers2(path):    # returns the outlier's index or -1 if not exist
    num_of_csbs = get_num_of_csbs_in_file(path)
    all_pq_tree= get_pqtree_from_file_and_index(path, 0)
    all_s_score = calculate_s_score(all_pq_tree)
    # #print("------------ all ", all_pq_tree_string, all_internal_nodes)

    outlier = -1
    num_outliers = 0
    for i in range(num_of_csbs):
        line_removed = remove_and_get_line(path, i)     # removing line at index i

        pq_tree = get_pqtree_from_file_and_index(path, 0)
        s_score = calculate_s_score(pq_tree)

        # #print("------------", i, pq_tree_string, internal_nodes)
        if s_score > all_s_score:
            outlier = i
            num_outliers += 1

        write_line_at_index(path, i, line_removed)     # writing the line at index i

    #print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! num of outliers:", num_outliers, "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n")

    # if num_outliers > 1:
    #     outlier = -1

    return outlier




def get_pqtrees_dict_with_outliers(path):
    outlier_index = find_outliers2(path)
    outlier_index = -1  # not allowing outliers

    if outlier_index == -1:     # there is no outliers
        num_of_csbs = get_num_of_csbs_in_file(path)
        pqtree_dict = {}
        for i in range(num_of_csbs):
            pqtree_dict[i] = {}
            pqtree_dict[i]["no_outlier"] = get_pqtree_from_file_and_index(path, i)
            pqtree_dict[i]["with_outlier"] = pqtree_dict[i]["no_outlier"]

        return pqtree_dict, outlier_index


    num_of_csbs = get_num_of_csbs_in_file(path)

    outlier_line_removed = remove_and_get_line(path, outlier_index)

    pqtree_dict = {}
    for i in range(num_of_csbs):
        pqtree_dict[i] = {}
        if i < outlier_index:
            pqtree_dict[i]["no_outlier"] = get_pqtree_from_file_and_index(path, i)
        if i > outlier_index:
            pqtree_dict[i]["no_outlier"] = get_pqtree_from_file_and_index(path, i-1)

    write_line_at_index(path, outlier_index, outlier_line_removed)  # writing the line at index i

    for i in range(num_of_csbs):
        pqtree_dict[i]["with_outlier"] = get_pqtree_from_file_and_index(path, i)

    return pqtree_dict, outlier_index



def get_pqtrees_dict_no_outliers(path):
    num_of_csbs = get_num_of_csbs_in_file(path)
    pqtree_dict = {}
    for i in range(num_of_csbs):
        pqtree_dict[i] = {}
        pqtree_dict[i]["no_outlier"] = get_pqtree_from_file_and_index(path, i)
        pqtree_dict[i]["with_outlier"] = pqtree_dict[i]["no_outlier"]

    return pqtree_dict, -1


def get_pqtrees_dict_no_outliers_tandem(path):
    num_of_csbs = get_num_of_csbs_in_file(path)
    pqtree_dict = {}
    for i in range(num_of_csbs):
        pqtree_dict[i] = {}
        pqtree_dict[i]["no_outlier"] = get_pqtree_from_file_and_index(path, i)
        pqtree_dict[i]["with_outlier"] = pqtree_dict[i]["no_outlier"]

    return pqtree_dict, -1




def add_penalty_for_outlier(MEM2, outlier_index, penalty):
    if outlier_index == -1:     # no outlier
        return MEM2

    for i in MEM2.keys():
        if i == outlier_index:
            MEM2[i] += penalty

    return MEM2


def get_tandem_dict(tandem_indeces):
    dict = {}
    for repeat in tandem_indeces:
        i = repeat[0]
        cog = repeat[1]
        if i not in dict:
            dict[i] = []
        dict[i].append(cog)
    return dict


def get_diff_repeats_cogs(cogs_first_csb,cogs_list2):
    diff = 0
    for cog in cogs_first_csb:
        if cog not in cogs_list2:
            diff += 1

    for cog in cogs_list2:
        if cog not in cogs_first_csb:
            diff += 1
    print("diffdiffdiffdiffdiffdiffdiffdiffdiff", cogs_first_csb, cogs_list2, diff)
    return diff


def get_diff_repeats_cogs2(cogs_first_csb,cogs_list2, tree_str, penaltyQ, penaltyP, path):
    diff = 0
    for cog in cogs_first_csb:
        if cog not in cogs_list2:
            penalty = get_penalty_tandem(cog, tree_str, penaltyQ, penaltyP, path)
            diff += penalty

    for cog in cogs_list2:
        if cog not in cogs_first_csb:
            penalty = get_penalty_tandem(cog, tree_str, penaltyQ, penaltyP, path)
            diff += penalty
    print("diffdiffdiffdiffdiffdiffdiffdiffdiff", cogs_first_csb, cogs_list2, diff)
    return diff


def add_penalty_for_tandem(MEM2, tandem_indeces, penalty):
    print(tandem_indeces)
    if len(tandem_indeces) == 0:
        return MEM2

    if tandem_indeces[0][0] != 0: #   the first csb has no tandem repeat
        for repeat in tandem_indeces:
            i = repeat[0]
            MEM2[i] += penalty
        return MEM2

    #   the first csb has tandem repeat
    tandem_dict = get_tandem_dict(tandem_indeces)
    for i in MEM2.keys():
        if i not in tandem_dict:
            MEM2[i] += penalty * get_diff_repeats_cogs(tandem_dict[0], [])
        else:
            MEM2[i] += penalty * get_diff_repeats_cogs(tandem_dict[0], tandem_dict[i])

    return MEM2

def get_penalty_tandem(cog, tree_str, penaltyQ, penaltyP, path):
    csbs_file_lines = open(path).readlines()
    first_csbs_list = csbs_file_lines[1].split('\t')[4].split(',')
    csbs_names_dict = {}
    for i in range(len(first_csbs_list)):
        csbs_names_dict[first_csbs_list[i][:-1]] = i

    cog_number = csbs_names_dict[cog]
    index = tree_str.index(str(cog_number))

    open_Q=0
    open_P=0
    for i in range(index,len(tree_str)):
        c = tree_str[i]
        if c == ')':
            open_P -= 1
            if open_P < 0:
                return penaltyP
        if c == ']':
            open_Q -= 1
            if open_Q < 0:
                return penaltyQ
        if c == '(':
            open_P += 1
        if c == '[':
            open_Q += 1


def add_penalty_for_tandem2(MEM2, tandem_indeces, tree_str, penaltyQ, penaltyP, path):
    print(tandem_indeces)
    if len(tandem_indeces) == 0:
        return MEM2

    if tandem_indeces[0][0] != 0: #   the first csb has no tandem repeat
        for repeat in tandem_indeces:
            i = repeat[0]
            cog = repeat[1]
            penalty = get_penalty_tandem(cog, tree_str, penaltyQ, penaltyP, path)
            MEM2[i] += penalty
        return MEM2

    #   the first csb has tandem repeat
    tandem_dict = get_tandem_dict(tandem_indeces)
    for i in MEM2.keys():
        if i not in tandem_dict:
            MEM2[i] += get_diff_repeats_cogs2(tandem_dict[0], [], tree_str, penaltyQ, penaltyP, path)
        else:
            MEM2[i] += get_diff_repeats_cogs2(tandem_dict[0], tandem_dict[i], tree_str, penaltyQ, penaltyP, path)

    return MEM2



def get_csbs_names_dict(csbs_file_path):
    csbs_file_lines = open(csbs_file_path).readlines()
    csbs_names_dict = {}
    for i in range(1,len(csbs_file_lines)):
        line = csbs_file_lines[i]
        csb = line.split('\t')[0]
        csbs_names_dict[i-1] = csb

    return csbs_names_dict


def pos_cogs_in_csb(csb):
    pos = 0
    for cog in csb:
        if cog[1] == 1:
            pos += 1
    return pos


def get_len_out_of_tree(pqtrees_dict):
    if "no_outlier" in pqtrees_dict:
        return span({pqtrees_dict["no_outlier"]})
    else:
        return span({pqtrees_dict["with_outlier"]})


def get_flipped_pattern(pattern):
    new_pattern = []
    pattern_len = len(pattern)
    for i in range(pattern_len):
        cog = pattern[pattern_len-1-i][0]
        sign = '+'
        if pattern[pattern_len-1-i][1] == '+':
            sign = '-'

        new_pattern.append((cog,sign))
    return new_pattern


# def check_tandem_and_create_temp_file(path):
#     csbs_file_lines = open(path).readlines()
#     first_line = csbs_file_lines[0]
#     fam_string = first_line[:]
#     tandem_indeces = []
#     index = 0
#     for line in csbs_file_lines[1:]:
#         csb = line.split('\t')[4].split(',')
#
#         cog_tandem = -1
#         for i in range(len(csb) - 1):
#             if csb[i] == csb[i + 1]:
#                 cog_tandem = csb[i]
#
#         if cog_tandem != -1:
#             line = line.replace(cog_tandem+","+cog_tandem, cog_tandem)
#             tandem_indeces.append(index)
#         index += 1
#
#         fam_string = fam_string + line
#
#     new_path = path.replace("new_families", "temp")
#     f = open(new_path, "a")
#     f.write(fam_string)
#     f.close()
#
#     return new_path, tandem_indeces



def check_tandem_and_create_temp_file(path):
    csbs_file_lines = open(path).readlines()
    first_line = csbs_file_lines[0]
    fam_string = first_line[:]
    tandem_indeces = []
    index = 0
    for line in csbs_file_lines[1:]:
        csb = line.split('\t')[4].split(',')

        for i in range(len(csb) - 1):
            if csb[i] == csb[i + 1]:
                cog_tandem = csb[i]
                line = line.replace(cog_tandem + "," + cog_tandem, cog_tandem)
                tandem_indeces.append((index,cog_tandem[:-1]))

        index += 1
        fam_string = fam_string + line

    new_path = path.replace("new_families", "temp")
    f = open(new_path, "a")
    f.write(fam_string)
    f.close()

    return new_path, tandem_indeces


def get_tandem_indices(path):
    csbs_file_lines = open(path).readlines()
    tandem_indeces = []
    index = 0
    for line in csbs_file_lines[1:]:
        csb = line.split('\t')[4].split(',')

        for i in range(len(csb) - 1):
            if csb[i] == csb[i + 1]:
                cog_tandem = csb[i]
                tandem_indeces.append((index,cog_tandem[:-1]))

        index += 1


    return tandem_indeces



def run_MEM4(path, d_T, d_S):
    path, tandem_indeces = check_tandem_and_create_temp_file(path)
    num_of_csbs = get_num_of_csbs_in_file(path)
    # pqtrees, outlier_index = get_pqtrees_dict_with_outliers(path)
    pqtrees, outlier_index = get_pqtrees_dict_no_outliers(path)



    pq_tree = 0
    MEM2 = {}
    # for i in range(num_of_csbs):
    i = 0
    # MEM2[i]={}
    pqtrees_list = [pqtrees[indx]["no_outlier"] for indx in range(num_of_csbs)]  # saved to use later for printing

    for j in range(num_of_csbs):
        # print("-----------------------------------CSBs: " + str(i) + " and   " + str(j))
        # #print("                                                                            csbs numbers: ", i, j)
        if i != outlier_index and j != outlier_index:
            pq_tree = deepcopy(pqtrees[i]["no_outlier"])
        else:
            pq_tree = deepcopy(pqtrees[i]["with_outlier"])


        pattern1 = get_pattern_by_file_and_index(path, j)
        A = main_algorithm(pq_tree, pattern1, d_T, d_S, {})
        # print(A)
        tree_span = get_len_out_of_tree(pqtrees[i])
        string_len = get_len_out_of_tree(pqtrees[j])
        # pos_values = [pos for pos in range(pos_cogs_in_csb(pattern))]
        mem1 = min([A[pq_tree][i][k_T][k_S][pos] for i in range(string_len-(tree_span-d_T)+1) for k_T in range(d_T+1) for k_S in range(d_S+1) for pos in range(tree_span-k_T+1)])


        pattern2 = get_flipped_pattern(pattern1)    # using the flipped csb as well
        A = main_algorithm(pq_tree, pattern2, d_T, d_S, {})
        # print(A)
        tree_span = get_len_out_of_tree(pqtrees[i])
        string_len = get_len_out_of_tree(pqtrees[j])
        # pos_values = [pos for pos in range(pos_cogs_in_csb(pattern))]
        mem2 = min([A[pq_tree][i][k_T][k_S][pos] for i in range(string_len-(tree_span-d_T)+1) for k_T in range(d_T+1) for k_S in range(d_S+1) for pos in range(tree_span-k_T+1)])


        MEM2[j] = min([mem1, mem2])
        # #TODO: check if the hole tree flipped
        # if pq_tree.type == "Q" and check_if_node_children_flip_penalties(pq_tree, 0, pq_tree.span - 1):
        #     MEM2[j] = MEM2[j] - qnode_flip_penalty*pq_tree.span

        # print(i, j, pq_tree, outlier_index, calculate_s_score(pq_tree), MEM2[i][j])
        # print(i, j, pq_tree, MEM2[j])

    MEM2 = add_penalty_for_outlier(MEM2, outlier_index, outlier_panalty)
    # MEM2 = add_penalty_for_tandem(MEM2, tandem_indeces, 0.5)
    MEM2 = add_penalty_for_tandem2(MEM2, tandem_indeces, str(pq_tree), 1.5, 1, path)

    return MEM2, pqtrees_list


def delete_all_temp_files():
    files = glob.glob(r'temp\*')
    for f in files:
        os.remove(f)


def get_all_MEM4_general_dicts(outlier_panal, bp_qnode_penal, qnode_flip_penal, d_T, d_S, delete_T_penal, delete_S_penal, jump_penal):
    global outlier_panalty
    global bp_qnode_penalty
    global qnode_flip_penalty
    global delete_T_penalty
    global delete_S_penalty
    global jump_penalty
    outlier_panalty = outlier_panal
    bp_qnode_penalty = bp_qnode_penal
    qnode_flip_penalty = qnode_flip_penal
    delete_T_penalty = delete_T_penal
    delete_S_penalty = delete_S_penal
    jump_penalty = jump_penal

    #print("\n\n---------- calculating MEM4 ---------------\n")
    csbs_files_paths = glob.glob(r'new_families\*.txt')

    MEM4_by_file_dict = {}
    pqtrees_list_dict = {}
    for path in csbs_files_paths:
        #print(path)
        file_name = path[13:].split(".")[0]
        print("\nGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG", file_name, "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n")
        MEM4_dict, pqtrees_list = run_MEM4(path, d_T, d_S)
        MEM4_by_file_dict[file_name] = MEM4_dict
        pqtrees_list_dict[file_name] = pqtrees_list

    delete_all_temp_files()
    return MEM4_by_file_dict, pqtrees_list_dict


# def print_all_MEM4_dicts(MEM4_by_file_dict):
#     for key in MEM4_by_file_dict:
#         #print("\n\n----------- ", key, " -------------")
#         for csb in MEM4_by_file_dict[key].keys():
#             #print(MEM4_by_file_dict[key][csb])


# def change_keys(MEM4_by_file_dict):
#     new_dict = {}
#     for path in MEM4_by_file_dict.keys():
#         new_dict[path] = {}
#         csbs_names_dict = get_csbs_names_dict("new_families/" + path + ".txt")
#         for key1 in MEM4_by_file_dict[path].keys():
#             new_dict[path][csbs_names_dict[key1]] = {}
#             for key2 in MEM4_by_file_dict[path].keys():
#                 # #print(csbs_names_dict)
#                 new_dict[path][csbs_names_dict[key1]][csbs_names_dict[key2]] = MEM4_by_file_dict[path][key1][key2]
#
#     return new_dict


def change_keys(MEM4_by_file_dict):
    new_dict = {}
    for path in MEM4_by_file_dict.keys():
        new_dict[path] = {}
        csbs_names_dict = get_csbs_names_dict("new_families/" + path + ".txt")
        for key1 in MEM4_by_file_dict[path].keys():
            new_dict[path][csbs_names_dict[key1]] = MEM4_by_file_dict[path][key1]

    return new_dict


def get_pqtrees_dict():
    csbs_files_paths = glob.glob(r'new_families\*.txt')

    pqtrees_dict = {}
    for path in csbs_files_paths:
        file_name = path[13:].split(".")[0]
        pqtrees, outlier_index = get_pqtrees_dict_with_outliers(path)
        pqtree = deepcopy(pqtrees[0]["no_outlier"])
        pqtrees_dict[file_name] = pqtree

    return pqtrees_dict


def main():
    # MEM4_by_file_dict = get_all_MEM4_dicts()
    # # print_all_MEM4_dicts(MEM4_by_file_dict)
    global outlier_panalty
    global bp_qnode_penalty
    global qnode_flip_penalty
    global delete_T_penalty
    global delete_S_penalty

    outlier_panalty = 3.6
    bp_qnode_penalty = 2.7
    qnode_flip_penalty = 1
    delete_T_penalty = 0
    delete_S_penalty = 0

    # print(run_MEM4("all_examples/onlyP.txt", 0, 0))

    # print(run_MEM4("all_examples/examplePQ_repeat.txt", 1, 1))
    path = "all_examples2/26238.txt"

    # all_pq_tree_string = get_pqtree_from_csbs_file_path_by_index(path, 0)
    # get_pqtree_from_file_and_index("all_examples/examplePQ_repeat.txt", 0)

    MEM4_dict = run_MEM4(path, 0, 0)
    # print(MEM4_dict)

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
