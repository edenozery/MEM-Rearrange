from pqtrees.pqtree_helpers.generate_s import IntervalHierarchy
from pqtrees.pqtree import PQTreeBuilder
from pqtrees.common_intervals.preprocess_find import common_k_indexed_with_singletons
from pqtrees.pqtree_helpers.reduce_intervals import ReduceIntervals
import glob
from itertools import chain, combinations


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


def get_pqtree_from_csbs_file_path_by_index(csbs_file_path, index):
    csbs_file_lines = open(csbs_file_path).readlines()
    first_csbs_list = csbs_file_lines[1].split('\t')[4].split(',')
    csbs_names_dict = {}
    for i in range(len(first_csbs_list)):
        csbs_names_dict[first_csbs_list[i][:-1]] = i

    perms = []
    for line in csbs_file_lines[1:]:
        perms.append(get_pattern_from_csb_line(line, csbs_names_dict))

    perms, index_dict, index_dict_reverse = get_perms_by_index(perms, index)

    common_intervals = common_k_indexed_with_singletons(*perms)

    ir_intervals = ReduceIntervals.reduce(common_intervals)
    s = IntervalHierarchy.from_irreducible_intervals(ir_intervals)

    pqtree = PQTreeBuilder._from_s(s)

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
    p_stack = [")"]
    q_stack = ["]"]
    children_list = []

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
        if len(children_list) == 2 and span != 2:
            return PQNode("P", children_list, span, direction, False, 0, 0, level)
        return PQNode("Q", children_list, span, direction, False, 0, 0, level)


def get_pattern_by_file_and_index(csbs_file_path, index):
    csbs_file_lines = open(csbs_file_path).readlines()
    first_csbs_list = csbs_file_lines[1].split('\t')[4].split(',')
    csbs_names_dict = {}
    for i in range(len(first_csbs_list)):
        csbs_names_dict[first_csbs_list[i][:-1]] = i

    line = csbs_file_lines[index+1]
    return get_pattern_list_from_csb_line(line, csbs_names_dict)


def get_pattern_directions_by_file_and_index(csbs_file_path, index):
    csbs_file_lines = open(csbs_file_path).readlines()
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


def powerset(iterable):
    "powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)"
    s = list(iterable)
    ps = list(chain.from_iterable(frozenset(combinations(s, r)) for r in range(len(s)+1)))
    ps = [frozenset(e) for e in ps]
    return ps


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
            penalties_sum += qnode_flip_penalty
        else:
            if child.type == "Q":
                penalties_sum += qnode_flip_penalty
    return penalties_sum


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

    return True


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


def p_mapping(pnode, i, e, csb, d_T, d_S, A):
    P = {}
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
                                if is_valid_entry(C, C_tag, N, k_T, k_S, y, pnode, i, e, sos_pos, csb, end_point):

                                    if len(C) == 1:
                                        y_start_index = end_point - length_L({y}, k_T, k_S) + 1
                                        P[C][C_tag][N][k_T][k_S][y][sos_pos] = A[y][y_start_index][k_T][k_S][sos_pos] + jump_violation_delta(C_tag,y)
                                    else:
                                        y_derivations_dict = get_y_derivations_dict(y, y_sign, end_point, k_T, k_S, sos_pos, A)

                                        case1 = get_P(P, C, C_tag, N, k_T, k_S-1, y, sos_pos) + delete_S_penalty
                                        case2 = float('inf')
                                        if len(y_derivations_dict.items()) > 0:
                                            case2 = min([get_P(P, C_minus_y, C_tag_minus_y, N_minus_y, k_T-deriv[0], k_S-deriv[1], z, sos_pos-deriv[2]) + score + dist_delta(N, y, z, pnode) + jump_violation_delta(C_tag,y) for (deriv, score) in y_derivations_dict.items() for z in C_minus_y])
                                        case3 = float('inf')
                                        P[C][C_tag][N][k_T][k_S][y][sos_pos] = min([case1, case2, case3])

                                else:
                                    P[C][C_tag][N][k_T][k_S][y][sos_pos] = float('inf')

    for k_T in range(d_T+1):
        for k_S in range(d_S+1):
            for pos in range(span({pnode}) + 1 - k_T):
                min_val = float('inf')
                for C in powerset(pnode.children):
                    if span({pnode}) - span(C) <= k_T:
                        min_val_C = min([get_P(P, C, C_tag, N, k_T, k_S, y, pos) + delete_T_penalty*span(set(pnode.children).difference(C)) for C_tag in powerset(C) for N in powerset(C) for y in C])
                        min_val = min([min_val, min_val_C])

                A[pnode][i][k_T][k_S][pos] = min_val
    return A


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

                A[qnode][i][k_T][k_S][pos] = min_val

    return A


def q_iteration(qnode, i, e, csb, d_T, d_S, Q, C, j, A, is_l):
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
                    N_r = {r: N.intersection(frozenset(C[:r])) for r in range(j)}
                    xj_sign = (xj not in N) * 2 - 1
                    end_point = end_point_E(C, i, k_T, k_S)
                    if is_valid_entry_Q(C, N, k_T, k_S, qnode, i, e, pos, xj, csb, end_point):
                        if j == 1:
                            xj_start_index = end_point - length_L({xj}, k_T, k_S) + 1
                            Q[j][k_T][k_S][N][pos] = A[xj][xj_start_index][k_T][k_S][pos]
                        else:
                            xj_derivations_dict = get_y_derivations_dict(xj, xj_sign, end_point, k_T, k_S, pos, A)

                            case1 = get_Q(Q, span_C[j], j, k_T, k_S-1, N, pos) + delete_S_penalty

                            case2 = float('inf')
                            if len(xj_derivations_dict.items()) > 0:
                                case2 = min([get_Q(Q, span_C[r], r, k_T - deriv[0] - span_r_to_j[r], k_S - deriv[1], N_r[r], pos - deriv[2]) + score + bp_qnode_penalty*dist_delta_Q2(N, xj, C[r-1], is_l) + delete_T_penalty*span_r_to_j[r] for r in [j-1] for (deriv, score) in xj_derivations_dict.items()]) # no deletions

                            case3 = float('inf')

                            Q[j][k_T][k_S][N][pos] = min([case1, case2, case3])

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
    flip_penalty = qnode_flip_penalty

    return flip_penalty - children_flipped_penalties


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
    return A


def iteration_internal_node(pqtree_node, csb, d_T, d_S, A):
    for i in range(len(csb)):
        e = end_point_E([pqtree_node], i, 0, d_S)

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
    else:
        for child in pqtree_node.children:
            A = main_algorithm(child, csb, d_T, d_S, A)
        A = iteration_internal_node(pqtree_node, csb, d_T, d_S, A)

    return A


def get_num_of_csbs_in_file(path):
    with open(path, 'r') as read_file:
        lines = read_file.readlines()
        num = len(lines) - 1
    return num


def get_pqtrees(path):
    num_of_csbs = get_num_of_csbs_in_file(path)
    pqtree_dict = {}
    for i in range(num_of_csbs):
        pqtree_dict[i] = get_pqtree_from_file_and_index(path, i)

    return pqtree_dict


def get_csbs_names_dict(csbs_file_path):
    csbs_file_lines = open(csbs_file_path).readlines()
    csbs_names_dict = {}
    for i in range(1,len(csbs_file_lines)):
        line = csbs_file_lines[i]
        csb = line.split('\t')[0]
        csbs_names_dict[i-1] = csb

    return csbs_names_dict


def get_len_out_of_tree(pqtree):
    return span({pqtree})


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


def run_MEM_first_CSB(path, d_T, d_S):  # runs all CSBs against the first CSB
    num_of_csbs = get_num_of_csbs_in_file(path)

    pqtrees = get_pqtrees(path)

    csbs_names_dict = get_csbs_names_dict(path)

    pq_tree = 0
    MEM2 = {}

    i = 0
    pqtrees_list = [pqtrees[indx] for indx in range(num_of_csbs)]  # saved to use later for printing
    print("PQ-tree:", pqtrees_list[0])
    print("MEM-Rearrange divergence scores:")

    for j in range(num_of_csbs):

        pq_tree = pqtrees[i]

        pattern1 = get_pattern_by_file_and_index(path, j)
        A = main_algorithm(pq_tree, pattern1, d_T, d_S, {})
        tree_span = get_len_out_of_tree(pq_tree)    # pqtrees[i]
        string_len = get_len_out_of_tree(pqtrees[j])
        mem1 = min([A[pq_tree][i][k_T][k_S][pos] for i in range(string_len-(tree_span-d_T)+1) for k_T in range(d_T+1) for k_S in range(d_S+1) for pos in range(tree_span-k_T+1)])

        pattern2 = get_flipped_pattern(pattern1)    # using the flipped csb as well
        A = main_algorithm(pq_tree, pattern2, d_T, d_S, {})
        tree_span = get_len_out_of_tree(pq_tree)    # pqtrees[i]
        string_len = get_len_out_of_tree(pqtrees[j])
        mem2 = min([A[pq_tree][i][k_T][k_S][pos] for i in range(string_len-(tree_span-d_T)+1) for k_T in range(d_T+1) for k_S in range(d_S+1) for pos in range(tree_span-k_T+1)])

        MEM2[csbs_names_dict[j]] = min([mem1, mem2])

    print(csbs_names_dict[i], ":     ", MEM2)

    return MEM2, pqtrees_list


def run_MEM_pairwise(path, d_T, d_S):  # pairwise: every CSB against every CSB
    num_of_csbs = get_num_of_csbs_in_file(path)

    pqtrees = get_pqtrees(path)

    csbs_names_dict = get_csbs_names_dict(path)

    MEM2 = {}

    pqtrees_list = [pqtrees[indx] for indx in range(num_of_csbs)]  # saved to use later for printing
    print("PQ-tree:", pqtrees_list[0])
    print("MEM-Rearrange divergence scores:")

    for i in range(num_of_csbs):
        MEM2[csbs_names_dict[i]] = {}
        for j in range(num_of_csbs):

            pq_tree = pqtrees[i]

            pattern1 = get_pattern_by_file_and_index(path, j)
            A = main_algorithm(pq_tree, pattern1, d_T, d_S, {})
            tree_span = get_len_out_of_tree(pq_tree)    # pqtrees[i]
            string_len = get_len_out_of_tree(pqtrees[j])
            mem1 = min([A[pq_tree][i][k_T][k_S][pos] for i in range(string_len-(tree_span-d_T)+1) for k_T in range(d_T+1) for k_S in range(d_S+1) for pos in range(tree_span-k_T+1)])

            pattern2 = get_flipped_pattern(pattern1)    # using the flipped csb as well
            A = main_algorithm(pq_tree, pattern2, d_T, d_S, {})
            tree_span = get_len_out_of_tree(pq_tree)    # pqtrees[i]
            string_len = get_len_out_of_tree(pqtrees[j])
            mem2 = min([A[pq_tree][i][k_T][k_S][pos] for i in range(string_len-(tree_span-d_T)+1) for k_T in range(d_T+1) for k_S in range(d_S+1) for pos in range(tree_span-k_T+1)])

            MEM2[csbs_names_dict[i]][csbs_names_dict[j]] = min([mem1, mem2])

        print(csbs_names_dict[i], ":     ", MEM2[csbs_names_dict[i]])

    return MEM2, pqtrees_list


def get_all_MEM_Rearrange_dicts(bp_qnode_penal, qnode_flip_penal, d_T, d_S, delete_T_penal, delete_S_penal, jump_penal):   # runs MEM-rearrange for each input file
    global bp_qnode_penalty
    global qnode_flip_penalty
    global delete_T_penalty
    global delete_S_penalty
    global jump_penalty
    bp_qnode_penalty = bp_qnode_penal
    qnode_flip_penalty = qnode_flip_penal
    delete_T_penalty = delete_T_penal
    delete_S_penalty = delete_S_penal
    jump_penalty = jump_penal

    print("\n---------- calculating MEM_general ---------------\n\n")

    csbs_files_paths = glob.glob(r'input_families\*.txt')

    MEM4_by_file_dict = {}
    pqtrees_list_dict = {}
    for path in csbs_files_paths:
        file_name = path[15:].split(".")[0]
        print("\n---------------------------------------", file_name, "---------------------------------------\n")

        # MEM4_dict, pqtrees_list = run_MEM_first_CSB(path, d_T, d_S)  # runs all CSBs against the first CSB
        MEM4_dict, pqtrees_list = run_MEM_pairwise(path, d_T, d_S)  # runs pairwise

        MEM4_by_file_dict[file_name] = MEM4_dict
        pqtrees_list_dict[file_name] = pqtrees_list

    return MEM4_by_file_dict






