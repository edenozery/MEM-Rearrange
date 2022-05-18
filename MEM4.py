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
    # #print(csbs_file_lines[1].split('\t')[4].split(','))
    first_csbs_list = csbs_file_lines[1].split('\t')[4].split(',')
    csbs_names_dict = {}
    for i in range(len(first_csbs_list)):
        csbs_names_dict[first_csbs_list[i][:-1]] = i

    perms = []
    for line in csbs_file_lines[1:]:
        perms.append(get_pattern_from_csb_line(line, csbs_names_dict))
    # #print(perms)

    perms, index_dict, index_dict_reverse = get_perms_by_index(perms, index)


    strs = {"".join(str(x) for x in p) for p in perms}

    common_intervals_trivial = trivial_common_k_with_singletons(*perms)
    common_intervals = common_k_indexed_with_singletons(*perms)

    ir_intervals = ReduceIntervals.reduce(common_intervals)
    s = IntervalHierarchy.from_irreducible_intervals(ir_intervals)

    pqtree = PQTreeBuilder._from_s(s)
    # #print(pqtree.to_parens())
    # #print(get_pqtree_by_dict(pqtree.to_parens(), index_dict_reverse))
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
    def __init__(self, type, children,span, direction, flipped, l, r):
        self.type = type
        self.children = children
        self.span = span
        self.direction = direction
        self.flipped = flipped
        self.l = l
        self.r = r

    def __str__(self):
        to_str = ""
        for child in self.children:
            to_str = to_str + str(child)

        if self.type == "P":
            to_str = "(" + to_str + ")"
        else:
            to_str = "[" + to_str + "]"
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
        return str(self.cog)


def get_pqtree_by_parens(parens, directions):
    # print("PPPPPPPPPPPPPPPPPPPPPPPPPPPP", parens)
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

    if parens[0]=="(":
        return PQNode("P", children_list, span, 1, False, 0, 0)
    else:
        return PQNode("Q", children_list, span, 1, False, 0, 0)


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
    print("pqtree: " + str(pq_tree))
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

























def iteration_MEM2_Leaf(leaf, csb, A):
    dict = {"+": 1, "-": -1}
    for i in range(len(csb) - 1, -1, -1):
        if csb[i][0] == leaf.cog:
            A[leaf][i] = 0
            leaf.l = i
            leaf.r = i
            if dict[csb[i][1]] != leaf.direction:
                leaf.flipped = True
                leaf.direction = -1
            else:
                leaf.direction = 1

        else:
            A[leaf][i] = float('inf')

    return A


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


def get_children_flipped_penalties(node):
    penalties_sum = 0
    for child in node.children:    # checking if all children were changed directions
        if type(child) is Leaf:
            penalties_sum += qnode_flip_penalty*1
        else:
            if child.type == "Q":
                penalties_sum += qnode_flip_penalty*child.span
    return penalties_sum


def check_if_node_children_flip_penalties(pnode, i, e):
    # index = i
    # flag1 = True
    # for j in range(len(pnode.children)): # checking children from left to right
    #     child = pnode.children[j]
    #     if child.l == index:
    #         index = child.r + 1
    #     else:
    #         flag1 = False
    #
    # if flag1 == False or index > e+1:
    #     flag1 = False

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


    # if flag1:
    #     children_dist = children_dist1
    #     children_order_string = children_order_string1
    # else:
    #     qnode.flipped = True
    #     children_dist = children_dist2
    #     children_order_string = children_order_string2

    all_children_flipped = check_if_all_children_flipped(pnode)
    all_children_changed_direction = check_if_all_children_changed_direction(pnode)

    if flipped and all_children_changed_direction and all_children_flipped:
        return True

    return False


def p_mapping(pnode, i, e, A):

    children_order_tree = [(str(i), pnode.children[i].direction) for i in range(len(pnode.children))]
    children_order_string = ""
    children_span_list = [child.span for child in pnode.children]
    index = i
    children_dist = 0
    while index <= e:
        next_index = -1
        for j in range(len(pnode.children)):
            child = pnode.children[j]
            if child.l == index:
                children_order_string = children_order_string + str(j)
                children_dist = children_dist + A[child][index]
                next_index = child.r + 1
        if next_index == -1 or next_index > e+1:
            return float('inf')
        index = next_index

    pnode.l = i
    pnode.r = e

    children_breakpoint = get_bp_from_children_orders(children_order_tree, children_order_string)
    jumping_penalty = get_jumping_penalty(children_order_tree, children_order_string, children_span_list)

    children_penalties = 0
    if check_if_node_children_flip_penalties(pnode, i, e):
        children_penalties = get_children_flipped_penalties(pnode)

    return children_breakpoint + children_dist + jumping_penalty - children_penalties


def q_mapping(qnode, i, e, A):
    children_order_tree = [(str(i), qnode.children[i].direction) for i in range(len(qnode.children))]
    children_order_string1 = ""
    children_order_string2 = ""

    index = i
    children_dist1 = 0
    flag1 = True
    for j in range(len(qnode.children)): # checking children from left to right
        child = qnode.children[j]
        if child.l == index:
            children_order_string1 = children_order_string1 + str(j)
            children_dist1 = children_dist1 + A[child][index]
            index = child.r + 1
        else:
            flag1 = False

    if flag1 == False or index > e+1:
        flag1 = False


    index = i
    children_dist2 = 0
    flag2 = True
    for j in range(len(qnode.children) - 1, -1, -1): # checking children from right to left
        child = qnode.children[j]
        if child.l == index:
            children_order_string2 = children_order_string2 + str(j)
            children_dist2 = children_dist2 + A[child][index]
            index = child.r + 1
        else:
            flag2 = False
    if flag2 == False or index > e+1:
        flag2 = False

    if flag1 == False and flag2 == False:
        return float('inf')


    children_dist = 0
    children_order_string = ""
    if flag1:
        children_dist = children_dist1
        children_order_string = children_order_string1
    else:
        qnode.flipped = True
        children_dist = children_dist2
        children_order_string = children_order_string2

    all_children_flipped = check_if_all_children_flipped(qnode)
    all_children_changed_direction = check_if_all_children_changed_direction(qnode)

    if all_children_changed_direction:
        qnode.direction = (-1)*qnode.direction

    qnode.l = i
    qnode.r = e

    flip_penalty = 0
    children_flipped_penalties = 0
    if qnode.flipped and all_children_changed_direction and all_children_flipped:
        flip_penalty = qnode_flip_penalty*qnode.span
        children_flipped_penalties = get_children_flipped_penalties(qnode)

    children_breakpoint = bp_qnode_penalty*get_bp_from_children_orders(children_order_tree, children_order_string)

    return children_dist + children_breakpoint + flip_penalty - children_flipped_penalties



def iteration_MEM2_internal_node(pqtree_node, csb, A):
    # print("CSBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB", csb)
    for i in range(len(csb) - 1, -1, -1):
        e = pqtree_node.span + i -1
        if e > len(csb)-1:
            continue
        if pqtree_node.type == "P":
            A[pqtree_node][i] = p_mapping(pqtree_node, i, e, A)
        else:
            # #print("QMAPPINGGGGGGGGGGG", pqtree_node, i)
            A[pqtree_node][i] = q_mapping(pqtree_node, i, e, A)
    return A


def calculate_MEM2(pqtree_node, csb, A):
    if type(pqtree_node) is Leaf:
        A[pqtree_node] = {}
        A = iteration_MEM2_Leaf(pqtree_node, csb, A)
    else:
        for child in pqtree_node.children:
            A[child] = {}
            A = calculate_MEM2(child, csb, A)
        A[pqtree_node] = {}
        A = iteration_MEM2_internal_node(pqtree_node, csb, A)

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

def add_penalty_for_outlier(MEM2, outlier_index, penalty):
    if outlier_index == -1:     # no outlier
        return MEM2

    for i in MEM2.keys():
        if i != outlier_index:
            MEM2[i][outlier_index] += penalty
            MEM2[outlier_index][i] += penalty

    return MEM2


def get_csbs_names_dict(csbs_file_path):
    csbs_file_lines = open(csbs_file_path).readlines()
    csbs_names_dict = {}
    for i in range(1,len(csbs_file_lines)):
        line = csbs_file_lines[i]
        csb = line.split('\t')[0]
        # #print("----------------------------------------------------------------------------------------------------", csb)
        csbs_names_dict[i-1] = csb

    return csbs_names_dict


def run_MEM4(path):
    num_of_csbs = get_num_of_csbs_in_file(path)
    pqtrees, outlier_index = get_pqtrees_dict_with_outliers(path)
    # #print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@", pqtrees[0]["no_outlier"])
    pq_tree = 0
    MEM2 = {}
    for i in range(num_of_csbs):
        MEM2[i]={}
        for j in range(num_of_csbs):
            # #print("                                                                            csbs numbers: ", i, j)
            if i != outlier_index and j != outlier_index:
                pq_tree = deepcopy(pqtrees[i]["no_outlier"])
            else:
                # print(f)
                pq_tree = deepcopy(pqtrees[i]["with_outlier"])

            pattern = get_pattern_by_file_and_index(path, j)
            A = calculate_MEM2(pq_tree, pattern, {})
            MEM2[i][j] = A[pq_tree][0]
            if pq_tree.type == "Q" and check_if_node_children_flip_penalties(pq_tree, 0, pq_tree.span - 1):
                MEM2[i][j] = MEM2[i][j] - qnode_flip_penalty*pq_tree.span
    MEM2 = add_penalty_for_outlier(MEM2, outlier_index, outlier_panalty)

    return MEM2


def get_all_MEM4_dicts(outlier_panal, bp_qnode_penal, qnode_flip_penal):
    global outlier_panalty
    global bp_qnode_penalty
    global qnode_flip_penalty
    outlier_panalty = outlier_panal
    bp_qnode_penalty = bp_qnode_penal
    qnode_flip_penalty = qnode_flip_penal

    #print("\n\n---------- calculating MEM4 ---------------\n")
    csbs_files_paths = glob.glob(r'all_examples2\*.txt')

    MEM4_by_file_dict = {}
    for path in csbs_files_paths:
        #print(path)
        MEM4_dict = run_MEM4(path)
        file_name = path[13:].split(".")[0]
        MEM4_by_file_dict[file_name] = MEM4_dict

    return MEM4_by_file_dict


# def print_all_MEM4_dicts(MEM4_by_file_dict):
#     for key in MEM4_by_file_dict:
#         #print("\n\n----------- ", key, " -------------")
#         for csb in MEM4_by_file_dict[key].keys():
#             #print(MEM4_by_file_dict[key][csb])


def change_keys(MEM4_by_file_dict):
    new_dict = {}
    for path in MEM4_by_file_dict.keys():
        new_dict[path] = {}
        csbs_names_dict = get_csbs_names_dict("new_families/" + path + ".txt")
        for key1 in MEM4_by_file_dict[path].keys():
            new_dict[path][csbs_names_dict[key1]] = {}
            for key2 in MEM4_by_file_dict[path].keys():
                # #print(csbs_names_dict)
                new_dict[path][csbs_names_dict[key1]][csbs_names_dict[key2]] = MEM4_by_file_dict[path][key1][key2]

    return new_dict



def main():
    # MEM4_by_file_dict = get_all_MEM4_dicts()
    # print_all_MEM4_dicts(MEM4_by_file_dict)
    global outlier_panalty
    global bp_qnode_penalty
    global qnode_flip_penalty
    outlier_panalty = 3.6
    bp_qnode_penalty = 2.7
    qnode_flip_penalty = 1
    # outlier_panalty = 0
    # bp_qnode_penalty = 0
    # qnode_flip_penalty = 0

    # print(run_MEM4("all_examples/onlyP.txt"))

    print(run_MEM4("all_examples2/26238.txt"))

# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
