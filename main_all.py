
from jaccard import *
from MEM4 import *
from breakpoint import *
from signed_breakpoint import get_all_signed_breakpoint_dicts
from MEM4_general import *
from MEM4_general_repeat import get_all_MEM4_general_repeat_dicts

from scipy.stats import pearsonr
from scipy.stats import spearmanr

import math
import random
import time
import numpy as np

import os

import copy


def print_correlations_from_dicts(dict1, dict2):
    sum_prsn_corr = 0
    sum_sprm_corr = 0
    dict1_keys_list = [key for key in dict1.keys()]
    dict2_keys_list = [key for key in dict2.keys()]
    for i in range(len(dict1_keys_list)):
        vector1 = []
        vector2 = []
        for j in range(len(dict1_keys_list)):
            vector1.append(1 - dict1[dict1_keys_list[i]][dict1_keys_list[j]])
            vector2.append(dict2[dict2_keys_list[i]][dict2_keys_list[j]])

        corr, _ = pearsonr(vector1, vector2)
        sum_prsn_corr = sum_prsn_corr + corr
        # #print('Pearsons correlation - breakpoint: %.3f' % corr)

        corr, _ = spearmanr(vector1, vector2)
        sum_sprm_corr = sum_sprm_corr + corr
        # #print('Spearmans correlation - breakpoint: %.3f' % corr)
        # #print('\n' + '-----------------------------------------------------------------------------------' + '\n')
    # #print('\n-----------------------------------------------------------------------------------')
    #print('AVG Pearsons correlations - MEM4: %.3f' % (sum_prsn_corr / len(dict1_keys_list)))
    #print('AVG Spearmans correlations - MEM4: %.3f' % (sum_sprm_corr / len(dict1_keys_list)))
    #print('-----------------------------------------------------------------------------------' + '\n')

    return (sum_prsn_corr / len(dict1_keys_list)), (sum_sprm_corr / len(dict1_keys_list))


def to_print(vector):
    s = ""
    for val in vector:
        s = s + "\t\t\t\t\t\t" + str(round(val, 3))
    return s



def print_correlations_from_dicts2(dict1, dict2, dict3):   # taking in account just the first csb distances

    dict1_keys_list = [key for key in dict1.keys()]
    dict2_keys_list = [key for key in dict2.keys()]
    # print("len of lists: ", len(dict1_keys_list), len(dict2_keys_list))
    i = 0
    vector_j,  vector_s, vector_n = [], [], []
    vector_m, vector_b = [], []
    for j in range(len(dict1_keys_list)):
        vector_j.append(dict1[dict1_keys_list[j]][0])   # Jaccard
        vector_s.append(dict1[dict1_keys_list[j]][1])   # Sorensen
        vector_n.append(dict1[dict1_keys_list[j]][2])   # Nestedness

        vector_m.append(dict2[dict2_keys_list[j]])   # MEM
        vector_b.append(dict3[dict2_keys_list[0]][dict2_keys_list[j]])   # SBP

    corr1_jm, _ = pearsonr(vector_j, vector_m)
    corr2_jm, _ = spearmanr(vector_j, vector_m)
    corr1_sm, _ = pearsonr(vector_s, vector_m)
    corr2_sm, _ = spearmanr(vector_s, vector_m)
    corr1_nm, _ = pearsonr(vector_n, vector_m)
    corr2_nm, _ = spearmanr(vector_n, vector_m)
    correlations_mem = [(corr1_jm, corr2_jm, corr1_sm), (corr2_sm, corr1_nm, corr2_nm)]

    corr1_jb, _ = pearsonr(vector_j, vector_b)
    corr2_jb, _ = spearmanr(vector_j, vector_b)
    corr1_sb, _ = pearsonr(vector_s, vector_b)
    corr2_sb, _ = spearmanr(vector_s, vector_b)
    corr1_nb, _ = pearsonr(vector_n, vector_b)
    corr2_nb, _ = spearmanr(vector_n, vector_b)
    correlations_sbp = [(corr1_jb, corr2_jb, corr1_sb), (corr2_sb, corr1_nb, corr2_nb)]

    #print('Pearsons correlations: %.3f' % corr1)
    #print('Spearmans correlations: %.3f' % corr2)
    #print('-----------------------------------------------------------------------------------' + '\n')
    print("MEM\t\t", to_print(vector_m))
    print("SBP\t\t", to_print(vector_b))
    print("Jaccard", to_print(vector_j))
    print("Sorensen", to_print(vector_s))
    print("Nestedness", to_print(vector_n))

    return correlations_mem, correlations_sbp



def get_correlations_from_dicts2(dict1, dict2, dict3):   # taking in account just the first csb distances

    dict1_keys_list = [key for key in dict1.keys()]
    dict2_keys_list = [key for key in dict2.keys()]
    # print("len of lists: ", len(dict1_keys_list), len(dict2_keys_list))
    i = 0
    vector_j,  vector_s, vector_n = [], [], []
    vector_m, vector_b = [], []
    for j in range(len(dict1_keys_list)):
        vector_j.append(dict1[dict1_keys_list[j]][0])   # Jaccard
        vector_s.append(dict1[dict1_keys_list[j]][1])   # Sorensen
        vector_n.append(dict1[dict1_keys_list[j]][2])   # Nestedness

        vector_m.append(dict2[dict2_keys_list[j]])   # MEM
        vector_b.append(dict3[dict2_keys_list[0]][dict2_keys_list[j]])   # SBP

    corr1_jm, _ = pearsonr(vector_j, vector_m)
    corr2_jm, _ = spearmanr(vector_j, vector_m)
    corr1_sm, _ = pearsonr(vector_s, vector_m)
    corr2_sm, _ = spearmanr(vector_s, vector_m)
    corr1_nm, _ = pearsonr(vector_n, vector_m)
    corr2_nm, _ = spearmanr(vector_n, vector_m)
    correlations_mem = [(corr1_jm, corr2_jm, corr1_sm), (corr2_sm, corr1_nm, corr2_nm)]

    corr1_jb, _ = pearsonr(vector_j, vector_b)
    corr2_jb, _ = spearmanr(vector_j, vector_b)
    corr1_sb, _ = pearsonr(vector_s, vector_b)
    corr2_sb, _ = spearmanr(vector_s, vector_b)
    corr1_nb, _ = pearsonr(vector_n, vector_b)
    corr2_nb, _ = spearmanr(vector_n, vector_b)
    correlations_sbp = [(corr1_jb, corr2_jb, corr1_sb), (corr2_sb, corr1_nb, corr2_nb)]


    return correlations_mem, correlations_sbp



def print_correlations_from_dicts2_edges(dict1, dict2, edges):   # taking in account just the first csb distances

    dict1_keys_list = [key for key in dict1.keys()]
    dict2_keys_list = [key for key in dict2.keys()]
    # print("len of lists: ", len(dict1_keys_list), len(dict2_keys_list))
    i = 0
    vector11,  vector12, vector13 = [], [], []
    vector2 = []
    # for j in range(len(dict1_keys_list)):
    #     vector11.append(dict1[dict1_keys_list[i]][dict1_keys_list[j]][0])
    #     vector12.append(dict1[dict1_keys_list[i]][dict1_keys_list[j]][1])
    #     vector13.append(dict1[dict1_keys_list[i]][dict1_keys_list[j]][2])
    #
    #     vector2.append(dict2[dict2_keys_list[i]][dict2_keys_list[j]])

    for edge in edges:
        vector11.append(dict1[edge[0]][edge[1]][0])
        vector12.append(dict1[edge[0]][edge[1]][1])
        vector13.append(dict1[edge[0]][edge[1]][2])

        vector2.append(dict2[edge[0]][edge[1]])

    corr11, _ = pearsonr(vector11, vector2)
    corr21, _ = spearmanr(vector11, vector2)

    corr12, _ = pearsonr(vector12, vector2)
    corr22, _ = spearmanr(vector12, vector2)

    corr13, _ = pearsonr(vector13, vector2)
    corr23, _ = spearmanr(vector13, vector2)

    #print('Pearsons correlations: %.3f' % corr1)
    #print('Spearmans correlations: %.3f' % corr2)
    #print('-----------------------------------------------------------------------------------' + '\n')

    return (corr11, corr12, corr13), (corr21, corr22, corr23)


def print_correlations_from_dicts3(dict1, dict2, edges):
    # #print(edges)
    sum_prsn_corr = 0
    sum_sprm_corr = 0
    dict1_keys_list = [key for key in dict1.keys()]
    dict2_keys_list = [key for key in dict2.keys()]
    num_of_correlations = 0
    for i in range(len(dict1_keys_list)):
        vector1 = []
        vector2 = []
        for edge in edges:
            if dict1_keys_list[i] == edge[0]:
                vector1.append(1 - dict1[edge[0]][edge[1]])
                vector2.append(dict2[edge[0]][edge[1]])

        if len(vector1) > 2:
            num_of_correlations += 1
            corr, _ = pearsonr(vector1, vector2)
            sum_prsn_corr = sum_prsn_corr + corr
            # #print('Pearsons correlation - breakpoint: %.3f' % corr)

            corr, _ = spearmanr(vector1, vector2)
            sum_sprm_corr = sum_sprm_corr + corr
            # #print('Spearmans correlation - breakpoint: %.3f' % corr)


    # #print('\n-----------------------------------------------------------------------------------')
    #print('AVG Pearsons correlations - MEM4: %.3f' % (sum_prsn_corr / num_of_correlations))
    #print('AVG Spearmans correlations - MEM4: %.3f' % (sum_sprm_corr / num_of_correlations))
    #print('-----------------------------------------------------------------------------------' + '\n')
    # #print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>", dict1.keys(), edges, num_of_correlations)
    return (sum_prsn_corr / num_of_correlations), (sum_sprm_corr / num_of_correlations)



def print_all_files_correlations_from_dicts_with_edges(files_dict1, files_dict2, all_edges_dict):
    print("\n----------------- correlatins --------------------")
    sum_prsn1, sum_sprm1 = 0, 0
    sum_prsn2, sum_sprm2 = 0, 0
    sum_prsn3, sum_sprm3 = 0, 0
    for key in files_dict1:
        print('\n\n---', key, '---')
        # print('\n\n---', key, '---', len(files_dict1[key]), len(files_dict2[key]))
        # prsn, sprm = print_correlations_from_dicts2(files_dict1[key], files_dict2[key])
        prsn, sprm = print_correlations_from_dicts2_edges(files_dict1[key], files_dict2[key], all_edges_dict[key])
        sum_prsn1 += prsn[0]
        sum_sprm1 += sprm[0]

        sum_prsn2 += prsn[1]
        sum_sprm2 += sprm[1]

        sum_prsn3 += prsn[2]
        sum_sprm3 += sprm[2]

        print("pearson:", prsn, "\nspearman:", sprm)


    print('-----------------------------------------------------------------------------------\n\n' + '\n')
    print('Jaccard - overall AVG Pearsons correlations: %.3f' % (sum_prsn1 / len(files_dict1.keys())))
    print('Jaccard - overall AVG Spearmans correlations: %.3f' % (sum_sprm1 / len(files_dict1.keys())))

    print('\nSorensen - overall AVG Pearsons correlations: %.3f' % (sum_prsn2 / len(files_dict1.keys())))
    print('Sorensen - overall AVG Spearmans correlations: %.3f' % (sum_sprm2 / len(files_dict1.keys())))

    print('\nNestedness - overall AVG Pearsons correlations: %.3f' % (sum_prsn3 / len(files_dict1.keys())))
    print('Nestedness - overall AVG Spearmans correlations: %.3f' % (sum_sprm3 / len(files_dict1.keys())))
    #print('-----------------------------------------------------------------------------------' + '\n')
    return (sum_prsn1 / len(files_dict1.keys()))


def print_pairs_trees(pqtrees_list, MEM_dict):
    MEM_keys = [key for key in MEM_dict.keys()]
    for i in range(len(pqtrees_list)):
        print(0, i, "  ", pqtrees_list[0],"  ", pqtrees_list[i], "  ", MEM_dict[MEM_keys[i]])


def print_cogs_info(file_name, cogs_info_dict):
    path = "new_families/" + file_name + ".txt"
    file_lines = open(path).readlines()
    first_csb = file_lines[1].split('\t')[4].split(',')
    num_instances = file_lines[1].split('\t')[3]
    print("\nfirst csb num of instances:", num_instances, "\n")
    # print("\nfirst csb:", first_csb, "\n")
    for i in range(len(first_csb)):
        print(i, first_csb[i][:-1], cogs_info_dict[first_csb[i][:-1]])


def get_cogs_info_dict():
    path = "all_examples/cog_info.txt"
    lines = open(path).readlines()
    cogs_info_dict = {}
    for line in lines:
        cogs_info_dict[line[:7]] = line[8:]

    return cogs_info_dict


def get_sorted_files_by_s_score(file_names):
    files_dict = {}

    for file_name in file_names:
        path = "new_families/" + file_name + ".txt"
        path, tandem_indeces = check_tandem_and_create_temp_file(path)
        pq_tree = get_pqtree_from_file_and_index(path, 0)
        s_score = calculate_s_score(pq_tree)
        files_dict[file_name] = s_score
        os.remove(path)

    results_sorted = sorted(files_dict.items(), key=lambda x:x[1], reverse=True)
    soted_file_names = [item[0] for item in results_sorted]

    return soted_file_names, files_dict


def calculate_q_score(pq_tree_str):
    return pq_tree_str.count('[')




def get_sorted_files_by_q_score(file_names):
    files_dict = {}

    for file_name in file_names:
        path = "new_families/" + file_name + ".txt"
        path, tandem_indeces = check_tandem_and_create_temp_file(path)
        pq_tree = get_pqtree_from_file_and_index(path, 0)
        q_score = calculate_q_score(str(pq_tree))
        files_dict[file_name] = q_score
        os.remove(path)

    results_sorted = sorted(files_dict.items(), key=lambda x:x[1], reverse=True)
    soted_file_names = [item[0] for item in results_sorted]

    return soted_file_names, files_dict



def print_corrs_by_q_score(files_q_score_corr_dict, files_q_score_corr_dict_sbp):
    for i in range(5):
        files_corrs = files_q_score_corr_dict[i]
        files_corrs_sbp = files_q_score_corr_dict_sbp[i]
        len_files_corrs = len(files_corrs)
        print("\n\n")
        print(len_files_corrs, "families with q-score: ", i)
        sum_prsn1_mem, sum_sprm1_mem = 0, 0
        sum_prsn2_mem, sum_sprm2_mem = 0, 0
        sum_prsn3_mem, sum_sprm3_mem = 0, 0
        sum_prsn1_sbp, sum_sprm1_sbp = 0, 0
        sum_prsn2_sbp, sum_sprm2_sbp = 0, 0
        sum_prsn3_sbp, sum_sprm3_sbp = 0, 0
        # for corrs_mem in files_corrs:
        for j in range(len(files_corrs)):
            corrs_mem = files_corrs[j]
            prsn_mem, sprm_mem = corrs_mem[0], corrs_mem[1]
            sum_prsn1_mem += prsn_mem[0]
            sum_sprm1_mem += sprm_mem[0]

            sum_prsn2_mem += prsn_mem[1]
            sum_sprm2_mem += sprm_mem[1]

            sum_prsn3_mem += prsn_mem[2]
            sum_sprm3_mem += sprm_mem[2]

            corrs_mem_sbp = files_corrs_sbp[j]
            prsn_sbp, sprm_sbp = corrs_mem_sbp[0], corrs_mem_sbp[1]
            sum_prsn1_sbp += prsn_sbp[0]
            sum_sprm1_sbp += sprm_sbp[0]

            sum_prsn2_sbp += prsn_sbp[1]
            sum_sprm2_sbp += sprm_sbp[1]

            sum_prsn3_sbp += prsn_sbp[2]
            sum_sprm3_sbp += sprm_sbp[2]

        print("\t\t\t\t\t\t\t\t MEM \t\t\t\t SBP \t\t\t diff")
        print("Jaccard - Pearsons correlations: ", round(sum_prsn1_mem / len_files_corrs, 3), "\t\t\t", round(sum_prsn1_sbp / len_files_corrs, 3), "\t\t\t", round((sum_prsn1_mem / len_files_corrs) - (sum_prsn1_sbp / len_files_corrs), 3))
        print("Jaccard - Spearmans correlations: ", round(sum_sprm1_mem / len_files_corrs, 3), "\t\t\t", round(sum_sprm1_sbp / len_files_corrs, 3), "\t\t\t", round((sum_sprm1_mem / len_files_corrs) - (sum_sprm1_sbp / len_files_corrs), 3))

        print("\nSorensen - Pearsons correlations: ", round(sum_prsn2_mem / len_files_corrs, 3), "\t\t\t", round(sum_prsn2_sbp / len_files_corrs, 3), "\t\t\t", round((sum_prsn2_mem / len_files_corrs) - (sum_prsn2_sbp / len_files_corrs), 3))
        print("Sorensen - Spearmans correlations: ", round(sum_sprm2_mem / len_files_corrs, 3), "\t\t\t", round(sum_sprm2_sbp / len_files_corrs, 3), "\t\t\t", round((sum_sprm2_mem / len_files_corrs) - (sum_sprm2_sbp / len_files_corrs), 3))

        print("\nNestedness - Pearsons correlations: ", round(sum_prsn3_mem / len_files_corrs, 3), "\t\t\t", round(sum_prsn3_sbp / len_files_corrs, 3), "\t\t\t", round((sum_prsn3_mem / len_files_corrs) - (sum_prsn3_sbp / len_files_corrs), 3))
        print("Nestedness - Spearmans correlations: ", round(sum_sprm3_mem / len_files_corrs, 3), "\t\t", round(sum_sprm3_sbp / len_files_corrs, 3), "\t\t\t", round((sum_sprm3_mem / len_files_corrs) - (sum_sprm3_sbp / len_files_corrs), 3))



def print_corrs_by_q_score2(files_q_score_corr_dict, files_q_score_corr_dict_sbp):
    sum_prsn1_mem, sum_sprm1_mem = 0, 0
    sum_prsn2_mem, sum_sprm2_mem = 0, 0
    sum_prsn3_mem, sum_sprm3_mem = 0, 0
    sum_prsn1_sbp, sum_sprm1_sbp = 0, 0
    sum_prsn2_sbp, sum_sprm2_sbp = 0, 0
    sum_prsn3_sbp, sum_sprm3_sbp = 0, 0

    files_corrs1 = files_q_score_corr_dict[1]
    files_corrs2 = files_q_score_corr_dict[2]
    files_corrs3 = files_q_score_corr_dict[3]
    files_corrs4 = files_q_score_corr_dict[4]
    len_files_corrs = len(files_corrs1) + len(files_corrs2) + len(files_corrs3) + len(files_corrs4)
    print("\n\n")
    print(len_files_corrs, "families")


    for i in [1,2,3,4]:
        files_corrs = files_q_score_corr_dict[i]
        files_corrs_sbp = files_q_score_corr_dict_sbp[i]

        # for corrs_mem in files_corrs:
        for j in range(len(files_corrs)):
            corrs_mem = files_corrs[j]
            prsn_mem, sprm_mem = corrs_mem[0], corrs_mem[1]
            sum_prsn1_mem += prsn_mem[0]
            sum_sprm1_mem += sprm_mem[0]

            sum_prsn2_mem += prsn_mem[1]
            sum_sprm2_mem += sprm_mem[1]

            sum_prsn3_mem += prsn_mem[2]
            sum_sprm3_mem += sprm_mem[2]

            corrs_mem_sbp = files_corrs_sbp[j]
            prsn_sbp, sprm_sbp = corrs_mem_sbp[0], corrs_mem_sbp[1]
            sum_prsn1_sbp += prsn_sbp[0]
            sum_sprm1_sbp += sprm_sbp[0]

            sum_prsn2_sbp += prsn_sbp[1]
            sum_sprm2_sbp += sprm_sbp[1]

            sum_prsn3_sbp += prsn_sbp[2]
            sum_sprm3_sbp += sprm_sbp[2]

    print("\t\t\t\t\t\t\t\t MEM \t\t\t\t SBP \t\t\t diff")
    print("Jaccard - Pearsons correlations: ", round(sum_prsn1_mem / len_files_corrs, 3), "\t\t\t", round(sum_prsn1_sbp / len_files_corrs, 3), "\t\t\t", round((sum_prsn1_mem / len_files_corrs) - (sum_prsn1_sbp / len_files_corrs), 3))
    print("Jaccard - Spearmans correlations: ", round(sum_sprm1_mem / len_files_corrs, 3), "\t\t\t", round(sum_sprm1_sbp / len_files_corrs, 3), "\t\t\t", round((sum_sprm1_mem / len_files_corrs) - (sum_sprm1_sbp / len_files_corrs), 3))

    print("\nSorensen - Pearsons correlations: ", round(sum_prsn2_mem / len_files_corrs, 3), "\t\t\t", round(sum_prsn2_sbp / len_files_corrs, 3), "\t\t\t", round((sum_prsn2_mem / len_files_corrs) - (sum_prsn2_sbp / len_files_corrs), 3))
    print("Sorensen - Spearmans correlations: ", round(sum_sprm2_mem / len_files_corrs, 3), "\t\t\t", round(sum_sprm2_sbp / len_files_corrs, 3), "\t\t\t", round((sum_sprm2_mem / len_files_corrs) - (sum_sprm2_sbp / len_files_corrs), 3))

    print("\nNestedness - Pearsons correlations: ", round(sum_prsn3_mem / len_files_corrs, 3), "\t\t\t", round(sum_prsn3_sbp / len_files_corrs, 3), "\t\t\t", round((sum_prsn3_mem / len_files_corrs) - (sum_prsn3_sbp / len_files_corrs), 3))
    print("Nestedness - Spearmans correlations: ", round(sum_sprm3_mem / len_files_corrs, 3), "\t\t", round(sum_sprm3_sbp / len_files_corrs, 3), "\t\t\t", round((sum_sprm3_mem / len_files_corrs) - (sum_sprm3_sbp / len_files_corrs), 3))



def print_corrs_by_q_score3(files_q_score_corr_dict, files_q_score_corr_dict_sbp):
    for i in range(5):
        files_corrs = files_q_score_corr_dict[i]
        files_corrs_sbp = files_q_score_corr_dict_sbp[i]
        len_files_corrs = len(files_corrs)
        print("\n\n")
        print(len_files_corrs, "families with q-score: ", i)
        sum_w = 0
        sum_prsn1_mem, sum_sprm1_mem = 0, 0
        sum_prsn2_mem, sum_sprm2_mem = 0, 0
        sum_prsn3_mem, sum_sprm3_mem = 0, 0
        sum_prsn1_sbp, sum_sprm1_sbp = 0, 0
        sum_prsn2_sbp, sum_sprm2_sbp = 0, 0
        sum_prsn3_sbp, sum_sprm3_sbp = 0, 0
        # for corrs_mem in files_corrs:
        for j in range(len(files_corrs)):
            sum_w
            corrs_mem = files_corrs[j]
            prsn_mem, sprm_mem = corrs_mem[0], corrs_mem[1]
            sum_prsn1_mem += prsn_mem[0]
            sum_sprm1_mem += sprm_mem[0]

            sum_prsn2_mem += prsn_mem[1]
            sum_sprm2_mem += sprm_mem[1]

            sum_prsn3_mem += prsn_mem[2]
            sum_sprm3_mem += sprm_mem[2]

            corrs_mem_sbp = files_corrs_sbp[j]
            prsn_sbp, sprm_sbp = corrs_mem_sbp[0], corrs_mem_sbp[1]
            sum_prsn1_sbp += prsn_sbp[0]
            sum_sprm1_sbp += sprm_sbp[0]

            sum_prsn2_sbp += prsn_sbp[1]
            sum_sprm2_sbp += sprm_sbp[1]

            sum_prsn3_sbp += prsn_sbp[2]
            sum_sprm3_sbp += sprm_sbp[2]

        print("\t\t\t\t\t\t\t\t MEM \t\t\t\t SBP \t\t\t diff")
        print("Jaccard - Pearsons correlations: ", round(sum_prsn1_mem / len_files_corrs, 3), "\t\t\t", round(sum_prsn1_sbp / len_files_corrs, 3), "\t\t\t", round((sum_prsn1_mem / len_files_corrs) - (sum_prsn1_sbp / len_files_corrs), 3))
        print("Jaccard - Spearmans correlations: ", round(sum_sprm1_mem / len_files_corrs, 3), "\t\t\t", round(sum_sprm1_sbp / len_files_corrs, 3), "\t\t\t", round((sum_sprm1_mem / len_files_corrs) - (sum_sprm1_sbp / len_files_corrs), 3))

        print("\nSorensen - Pearsons correlations: ", round(sum_prsn2_mem / len_files_corrs, 3), "\t\t\t", round(sum_prsn2_sbp / len_files_corrs, 3), "\t\t\t", round((sum_prsn2_mem / len_files_corrs) - (sum_prsn2_sbp / len_files_corrs), 3))
        print("Sorensen - Spearmans correlations: ", round(sum_sprm2_mem / len_files_corrs, 3), "\t\t\t", round(sum_sprm2_sbp / len_files_corrs, 3), "\t\t\t", round((sum_sprm2_mem / len_files_corrs) - (sum_sprm2_sbp / len_files_corrs), 3))

        print("\nNestedness - Pearsons correlations: ", round(sum_prsn3_mem / len_files_corrs, 3), "\t\t\t", round(sum_prsn3_sbp / len_files_corrs, 3), "\t\t\t", round((sum_prsn3_mem / len_files_corrs) - (sum_prsn3_sbp / len_files_corrs), 3))
        print("Nestedness - Spearmans correlations: ", round(sum_sprm3_mem / len_files_corrs, 3), "\t\t", round(sum_sprm3_sbp / len_files_corrs, 3), "\t\t\t", round((sum_sprm3_mem / len_files_corrs) - (sum_sprm3_sbp / len_files_corrs), 3))




def get_sorted_files_by_corr(files_dict1, files_dict2, files_dict3):
    files_dict = {}

    for key in files_dict1:
        correlations_mem, correlations_sbp = get_correlations_from_dicts2(files_dict1[key], files_dict2[key], files_dict3[key])
        prsn_mem, sprm_mem = correlations_mem[0], correlations_mem[1]
        score = prsn_mem[0]
        files_dict[key] = score

    results_sorted = sorted(files_dict.items(), key=lambda x:x[1], reverse=True)
    soted_file_names = [item[0] for item in results_sorted]

    return soted_file_names


def print_tandem_list(key):
    path = "new_families/" + key + ".txt"
    tandem_indeces = get_tandem_indices(path)
    print("\ntandem repeats: ", tandem_indeces,"\n")



def print_all_files_correlations_from_dicts(files_dict1, files_dict2, files_dict3, pqtrees_dict):  # 1 - Jaccard  ;  2 - MEM  ;  3 - SBP
    print("\n----------------- correlatins --------------------")
    sum_prsn1_mem, sum_sprm1_mem = 0, 0
    sum_prsn2_mem, sum_sprm2_mem = 0, 0
    sum_prsn3_mem, sum_sprm3_mem = 0, 0

    sum_prsn1_sbp, sum_sprm1_sbp = 0, 0
    sum_prsn2_sbp, sum_sprm2_sbp = 0, 0
    sum_prsn3_sbp, sum_sprm3_sbp = 0, 0

    cogs_info_dict = get_cogs_info_dict()

    files_sorted_by_s_score, s_score_dict = get_sorted_files_by_s_score([k for k in files_dict1.keys()])

    files_sorted_by_q_score, q_score_dict = get_sorted_files_by_q_score([k for k in files_dict1.keys()])

    files_sorted_by_corr = get_sorted_files_by_corr(files_dict1, files_dict2, files_dict3)


    files_q_score_corr_dict = {0: [], 1: [], 2: [], 3: [], 4: []}
    files_q_score_corr_dict_sbp = {0: [], 1: [], 2: [], 3: [], 4: []}

    sum_w = 0
    # for key in files_dict1:
    for key in files_sorted_by_corr:
        print('\n\n---', key, '---')
        print("s-score", s_score_dict[key])
        print("Q-score", q_score_dict[key])
        print_cogs_info(key, cogs_info_dict)
        # print(pqtrees_dict[key])
        print_pairs_trees(pqtrees_dict[key], files_dict2[key])
        print_tandem_list(key)
        len_family = len(pqtrees_dict[key])
        sum_w += len_family
        # print('\n\n---', key, '---', len(files_dict1[key]), len(files_dict2[key]))
        correlations_mem, correlations_sbp = print_correlations_from_dicts2(files_dict1[key], files_dict2[key], files_dict3[key])
        prsn_mem, sprm_mem = correlations_mem[0], correlations_mem[1]
        sum_prsn1_mem += prsn_mem[0] * len_family
        sum_sprm1_mem += sprm_mem[0] * len_family

        sum_prsn2_mem += prsn_mem[1] * len_family
        sum_sprm2_mem += sprm_mem[1] * len_family

        sum_prsn3_mem += prsn_mem[2] * len_family
        sum_sprm3_mem += sprm_mem[2] * len_family

        print("\nMEM - pearson:", prsn_mem, "\nspearman:", sprm_mem)

        files_q_score_corr_dict[q_score_dict[key]].append((prsn_mem, sprm_mem))


        prsn_sbp, sprm_sbp = correlations_sbp[0], correlations_sbp[1]
        sum_prsn1_sbp += prsn_sbp[0] * len_family
        sum_sprm1_sbp += sprm_sbp[0] * len_family

        sum_prsn2_sbp += prsn_sbp[1] * len_family
        sum_sprm2_sbp += sprm_sbp[1] * len_family

        sum_prsn3_sbp += prsn_sbp[2] * len_family
        sum_sprm3_sbp += sprm_sbp[2] * len_family

        print("\nSBP - pearson:", prsn_sbp, "\nspearman:", sprm_sbp)

        files_q_score_corr_dict_sbp[q_score_dict[key]].append((prsn_sbp, sprm_sbp))

    # MEM
    print('\n\n\n------------------------------------------MEM-----------------------------------------\n')
    print('Jaccard - overall AVG Pearsons correlations: %.3f' % (sum_prsn1_mem / sum_w))
    print('Jaccard - overall AVG Spearmans correlations: %.3f' % (sum_sprm1_mem / sum_w))

    print('\nSorensen - overall AVG Pearsons correlations: %.3f' % (sum_prsn2_mem / sum_w))
    print('Sorensen - overall AVG Spearmans correlations: %.3f' % (sum_sprm2_mem / sum_w))

    print('\nNestedness - overall AVG Pearsons correlations: %.3f' % (sum_prsn3_mem / sum_w))
    print('Nestedness - overall AVG Spearmans correlations: %.3f' % (sum_sprm3_mem / sum_w))

    #   SBP
    print('\n\n------------------------------------------SBP-----------------------------------------\n')
    print('Jaccard - overall AVG Pearsons correlations: %.3f' % (sum_prsn1_sbp / sum_w))
    print('Jaccard - overall AVG Spearmans correlations: %.3f' % (sum_sprm1_sbp / sum_w))

    print('\nSorensen - overall AVG Pearsons correlations: %.3f' % (sum_prsn2_sbp / sum_w))
    print('Sorensen - overall AVG Spearmans correlations: %.3f' % (sum_sprm2_sbp / sum_w))

    print('\nNestedness - overall AVG Pearsons correlations: %.3f' % (sum_prsn3_sbp / sum_w))
    print('Nestedness - overall AVG Spearmans correlations: %.3f' % (sum_sprm3_sbp / sum_w))

    print("\n\n")

    print_corrs_by_q_score(files_q_score_corr_dict, files_q_score_corr_dict_sbp)
    # print_corrs_by_q_score2(files_q_score_corr_dict, files_q_score_corr_dict_sbp)
    print_corrs_by_q_score3(files_q_score_corr_dict, files_q_score_corr_dict_sbp)

    #print('-----------------------------------------------------------------------------------' + '\n')
    # return (sum_prsn1 / len(files_dict1.keys()))


def print_all_files_correlations_from_dicts3(files_dict1, files_dict2, files_dict3, pqtrees_dict):  # 1 - Jaccard  ;  2 - MEM  ;  3 - SBP
    print("\n----------------- correlatins --------------------")
    sum_prsn1_mem, sum_sprm1_mem = 0, 0
    sum_prsn2_mem, sum_sprm2_mem = 0, 0
    sum_prsn3_mem, sum_sprm3_mem = 0, 0

    sum_prsn1_sbp, sum_sprm1_sbp = 0, 0
    sum_prsn2_sbp, sum_sprm2_sbp = 0, 0
    sum_prsn3_sbp, sum_sprm3_sbp = 0, 0

    for key in files_dict1:
        # print('\n\n---', key, '---')
        # print(pqtrees_dict[key])
        print_pairs_trees(pqtrees_dict[key], files_dict2[key])
        # print('\n\n---', key, '---', len(files_dict1[key]), len(files_dict2[key]))
        correlations_mem, correlations_sbp = print_correlations_from_dicts2(files_dict1[key], files_dict2[key], files_dict3[key])
        prsn_mem, sprm_mem = correlations_mem[0], correlations_mem[1]
        sum_prsn1_mem += prsn_mem[0]
        sum_sprm1_mem += sprm_mem[0]

        sum_prsn2_mem += prsn_mem[1]
        sum_sprm2_mem += sprm_mem[1]

        sum_prsn3_mem += prsn_mem[2]
        sum_sprm3_mem += sprm_mem[2]

        # print("\nMEM - pearson:", prsn_mem, "\nspearman:", sprm_mem)


        prsn_sbp, sprm_sbp = correlations_sbp[0], correlations_sbp[1]
        sum_prsn1_sbp += prsn_sbp[0]
        sum_sprm1_sbp += sprm_sbp[0]

        sum_prsn2_sbp += prsn_sbp[1]
        sum_sprm2_sbp += sprm_sbp[1]

        sum_prsn3_sbp += prsn_sbp[2]
        sum_sprm3_sbp += sprm_sbp[2]

        # print("\nSBP - pearson:", prsn_sbp, "\nspearman:", sprm_sbp)



    # MEM
    print('\n\n\n------------------------------------------MEM-----------------------------------------\n')
    print('Jaccard - overall AVG Pearsons correlations: %.3f' % (sum_prsn1_mem / len(files_dict1.keys())))
    print('Jaccard - overall AVG Spearmans correlations: %.3f' % (sum_sprm1_mem / len(files_dict1.keys())))

    print('\nSorensen - overall AVG Pearsons correlations: %.3f' % (sum_prsn2_mem / len(files_dict1.keys())))
    print('Sorensen - overall AVG Spearmans correlations: %.3f' % (sum_sprm2_mem / len(files_dict1.keys())))

    print('\nNestedness - overall AVG Pearsons correlations: %.3f' % (sum_prsn3_mem / len(files_dict1.keys())))
    print('Nestedness - overall AVG Spearmans correlations: %.3f' % (sum_sprm3_mem / len(files_dict1.keys())))

    #   SBP
    print('\n\n------------------------------------------SBP-----------------------------------------\n')
    print('Jaccard - overall AVG Pearsons correlations: %.3f' % (sum_prsn1_sbp / len(files_dict1.keys())))
    print('Jaccard - overall AVG Spearmans correlations: %.3f' % (sum_sprm1_sbp / len(files_dict1.keys())))

    print('\nSorensen - overall AVG Pearsons correlations: %.3f' % (sum_prsn2_sbp / len(files_dict1.keys())))
    print('Sorensen - overall AVG Spearmans correlations: %.3f' % (sum_sprm2_sbp / len(files_dict1.keys())))

    print('\nNestedness - overall AVG Pearsons correlations: %.3f' % (sum_prsn3_sbp / len(files_dict1.keys())))
    print('Nestedness - overall AVG Spearmans correlations: %.3f' % (sum_sprm3_sbp / len(files_dict1.keys())))

    #print('-----------------------------------------------------------------------------------' + '\n')

    return [((sum_prsn1_mem / len(files_dict1.keys())),(sum_sprm1_mem / len(files_dict1.keys()))), ((sum_prsn2_mem / len(files_dict1.keys())), (sum_sprm2_mem / len(files_dict1.keys()))), ((sum_prsn3_mem / len(files_dict1.keys())), (sum_sprm3_mem / len(files_dict1.keys())))]


def print_all_files_correlations_from_dicts_by_edges(files_dict1, files_dict2, files_edges_dict):
    # #print("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1", files_dict1.keys())
    sum_prsn = 0
    sum_sprm = 0
    num_edges = 0
    for key in files_dict1:
        #print('\n------------------------------------------', key, '-----------------------------------------')
        prsn, sprm = print_correlations_from_dicts3(files_dict1[key], files_dict2[key], files_edges_dict[key])
        sum_prsn += prsn
        sum_sprm += sprm
        num_edges += len(files_edges_dict[key])

    #print('-----------------------------------------------------------------------------------\n\n' + '\n')
    #print('overall AVG Pearsons correlations: %.3f' % (sum_prsn / num_edges))
    #print('overall AVG Spearmans correlations: %.3f' % (sum_sprm / num_edges))
    #print('-----------------------------------------------------------------------------------' + '\n')


def get_minimum_spanning_tree_jaccard(edges, jaccard_index_dict):
    new_edges = []
    for node in jaccard_index_dict.keys():
        min_weight = float('inf')
        minimum_in_edges = []
        for edge in edges:
            if edge[1] == node:
                weight = (1 - jaccard_index_dict[edge[0]][edge[1]])
                if weight == min_weight:
                    minimum_in_edges.append(edge)
                if weight < min_weight:
                    min_weight = weight
                    minimum_in_edges = [edge]

        new_edges.extend(minimum_in_edges)

    return new_edges


def get_minimum_spanning_tree(edges, dict_bp_MEM):
    new_edges = []
    for node in dict_bp_MEM.keys():
        min_weight = float('inf')
        minimum_in_edges = []
        for edge in edges:
            if edge[1] == node:
                weight = dict_bp_MEM[edge[0]][edge[1]]
                if weight == min_weight:
                    minimum_in_edges.append(edge)
                if weight < min_weight:
                    min_weight = weight
                    minimum_in_edges = [edge]

        new_edges.extend(minimum_in_edges)

    return new_edges


def calculate_edges_diff(jaccard_index_dict, bp_MEM_dict, edges):
    # jaccard_index_dict, nodes, edges = get_jaccard_index_nodes_edges_by_path(path_name + "_instances.fasta")
    # breakpoint_dict = get_breakpoint_from_file(path_name + ".txt")

    new_edges_jaccard = get_minimum_spanning_tree_jaccard(edges, jaccard_index_dict)
    new_edges_breakpoint = get_minimum_spanning_tree(edges, bp_MEM_dict)

    difference_edges1 = set(new_edges_breakpoint).difference(set(new_edges_jaccard))
    difference_edges2 = set(new_edges_jaccard).difference(set(new_edges_breakpoint))
    return len(difference_edges1) + len(difference_edges2)


def get_in_edges(node, edges):
    in_edges = []
    for edge in edges:
        if edge[1] == node:
            in_edges.append(edge)
    return in_edges


def get_diff_by_edges(nodes, new_edges_jaccard, new_edges_breakpoint, net_edges):
    total_diff = 0
    for node in nodes:
        in_edges_net = get_in_edges(node, net_edges)
        if len(in_edges_net) > 0:
            in_edges_jaccard = get_in_edges(node, new_edges_jaccard)
            in_edges_bp = get_in_edges(node, new_edges_breakpoint)
            # similarity = len(set(in_edges_jaccard).intersection(set(in_edges_bp))) / len(in_edges_jaccard)
            # difference = len(in_edges_bp) / len(net_edges)
            similarity = len(set(in_edges_jaccard).intersection(set(in_edges_bp)))
            difference = len(set(in_edges_bp).difference(set(in_edges_jaccard)))
            total_diff += similarity - difference
            # if len(set(in_edges_jaccard).difference(set(in_edges_bp))) == 0:
            #     total_diff += math.log2(len(in_edges_net) / len(in_edges_bp))
            # else:
            #     total_diff += math.log2(1 - (len(in_edges_bp) / len(in_edges_net)))

    return total_diff


def calculate_edges_diff2(jaccard_index_dict, bp_MEM_dict, edges):

    new_edges_jaccard = get_minimum_spanning_tree_jaccard(edges, jaccard_index_dict)
    new_edges_breakpoint = get_minimum_spanning_tree(edges, bp_MEM_dict)
    # new_edges_breakpoint = new_edges_jaccard    # for getting jaccard vs jaccard

    diff = get_diff_by_edges(jaccard_index_dict.keys(), new_edges_jaccard, new_edges_breakpoint, edges)

    # difference_edges1 = set(new_edges_breakpoint).difference(set(new_edges_jaccard))
    # difference_edges2 = set(new_edges_jaccard).difference(set(new_edges_breakpoint))
    return diff


def get_trees_difference(jacard_index_by_file_dict, bp_MEM_by_file_dict, all_edges_dict):
    files_paths = glob.glob(r'all_examples\*.txt')
    edges_diff = 0
    for path in files_paths:
        path_name = path[13:].split(".")[0]
        edges_diff += calculate_edges_diff2(jacard_index_by_file_dict[path_name], bp_MEM_by_file_dict[path_name], all_edges_dict[path_name])
    return edges_diff, (edges_diff/len(files_paths))


def get_trees_difference_by_files(jacard_index_by_file_dict, bp_MEM_by_file_dict, all_edges_dict, files_names):
    edges_diff = 0
    for path_name in files_names:
        edges_diff += calculate_edges_diff2(jacard_index_by_file_dict[path_name], bp_MEM_by_file_dict[path_name], all_edges_dict[path_name])
    return edges_diff



def run_grid_search(jacard_index_by_file_dict, all_edges_dict):
    outlier_panalties = [0]
    bp_qnode_penalties = [1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3, 3.25, 3.5, 3.75, 4]
    flip_penalties = [0]
    output_dict = {}

    for penalty1 in outlier_panalties:
        for penalty2 in bp_qnode_penalties:
            for penalty3 in flip_penalties:
                # MEM4_by_file_dict = get_all_MEM4_dicts(penalty1, penalty2,0)
                MEM4_by_file_dict = get_all_MEM4_general_dicts(penalty1, penalty2, penalty3, 0, 0, 0, 0)
                MEM4_by_file_dict = change_keys(MEM4_by_file_dict)
                corr = print_all_files_correlations_from_dicts(jacard_index_by_file_dict, MEM4_by_file_dict)
                # trees_diff, trees_diff_AVG = get_trees_difference(jacard_index_by_file_dict, MEM4_by_file_dict, all_edges_dict)
                output_dict[(penalty1, penalty2)] = corr

    print("\n\n")
    print(output_dict)

def argmax(lst):
    index = 0
    max_val = 0
    for i in range(len(lst)):
        if lst[i] > max_val:
            index = i
            max_val = lst[i]
    return index, max_val



def print_sorted_results(corr_dict, corr_name):
    # results_sorted = sorted(corr_dict.items(), key=lambda x:x[1])
    # top_20 = [results_sorted[i][0] for i in range(-1, -100, -1)]
    # print(corr_name, "\n", top_20, "\n\n")

    results_sorted = sorted(corr_dict.items(), key=lambda x:x[1], reverse=True)
    print(corr_name, "\n", results_sorted, "\n\n")



def run_grid_search3(jacard_index_by_file_dict, signed_breakpoint_by_file_dict):
    outlier_panalties = [0]
    # bp_qnode_penalties = [0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2]
    bp_qnode_penalties = [1.5]
    # flip_penalties = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
    flip_penalties = [0.5]
    jump_penalties = [0, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2, 2.2, 2.4, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 3.8, 4]

    iter = 1
    len_penals_comb = len(outlier_panalties) * len(bp_qnode_penalties) * len(flip_penalties) * len(jump_penalties)
    jaccard_pearson, jaccard_spearman, sorensen_pearson, sorensen_spearman, nestedness_pearson, nestedness_spearman = {}, {}, {}, {}, {}, {}
    for penalty1 in outlier_panalties:
        for penalty2 in bp_qnode_penalties:
            for penalty3 in flip_penalties:
                for penalty4 in jump_penalties:
                    print(iter, "/", len_penals_comb)
                    iter += 1

                    MEM4_general_by_file_dict, pqtrees_dict = get_all_MEM4_general_dicts(penalty1, penalty2, penalty3, 0, 0, 0, 0, penalty4)
                    MEM4_general_by_file_dict = change_keys(MEM4_general_by_file_dict)

                    results = print_all_files_correlations_from_dicts3(jacard_index_by_file_dict, MEM4_general_by_file_dict,
                                                            signed_breakpoint_by_file_dict, pqtrees_dict)

                    jaccard_pearson[(penalty1, penalty2, penalty3, penalty4)] = results[0][0]
                    jaccard_spearman[(penalty1, penalty2, penalty3, penalty4)] = results[0][1]

                    sorensen_pearson[(penalty1, penalty2, penalty3, penalty4)] = results[1][0]
                    sorensen_spearman[(penalty1, penalty2, penalty3, penalty4)] = results[1][1]

                    nestedness_pearson[(penalty1, penalty2, penalty3, penalty4)] = results[2][0]
                    nestedness_spearman[(penalty1, penalty2, penalty3, penalty4)] = results[2][1]


                    # trees_diff, trees_diff_AVG = get_trees_difference(jacard_index_by_file_dict, MEM4_by_file_dict, all_edges_dict)
                    # output_dict[(penalty1, penalty2)] = corr

    print("\n\n")

    print_sorted_results(jaccard_pearson, "jaccard_pearson")
    print_sorted_results(jaccard_spearman, "jaccard_spearman")
    print_sorted_results(sorensen_pearson, "sorensen_pearson")
    print_sorted_results(sorensen_spearman, "sorensen_spearman")
    print_sorted_results(nestedness_pearson, "nestedness_pearson")
    print_sorted_results(nestedness_spearman, "nestedness_spearman")





def print_all_files_correlations_from_dicts2(jacard_index_by_file_dict, MEM4_by_file_dict, files_names_samples):
    sum_prsn = 0
    for key in files_names_samples:
        print('\n------------------------------------------', key, '-----------------------------------------')
        prsn, sprm = print_correlations_from_dicts2(jacard_index_by_file_dict[key], MEM4_by_file_dict[key])
        sum_prsn += prsn

    return (sum_prsn / len(files_names_samples))




def get_top_20_results(jacard_index_by_file_dict, all_penalties_MEM4_dict, all_edges_dict):
    files_names_samples = random.sample(jacard_index_by_file_dict.keys(), 10)
    examples_names = [key for key in jacard_index_by_file_dict.keys()]
    test_set = set(examples_names).difference(set(files_names_samples))
    results_dict = {}
    for penalty1, penalty2, penalty3 in all_penalties_MEM4_dict.keys():
        corr = print_all_files_correlations_from_dicts2(jacard_index_by_file_dict, all_penalties_MEM4_dict[(penalty1, penalty2, penalty3)], files_names_samples)
        # corr = get_trees_difference_by_files(jacard_index_by_file_dict, all_penalties_MEM4_dict[(penalty1, penalty2)], all_edges_dict, files_names_samples)
        results_dict[(penalty1, penalty2, penalty3)] = corr

    results_sorted = sorted(results_dict.items(), key=lambda x:x[1])
    top_20 = [results_sorted[i][0] for i in range(-1, -100, -1)]
    # top_20 = [results_sorted[i][0] for i in range(-1, -5, -1)]

    # check on 10 test examples
    top_20_reults_dict = {}
    for penalty1, penalty2, penalty3 in top_20:
        corr = print_all_files_correlations_from_dicts2(jacard_index_by_file_dict, all_penalties_MEM4_dict[(penalty1, penalty2, penalty3)], files_names_samples)
        # corr = get_trees_difference_by_files(jacard_index_by_file_dict, all_penalties_MEM4_dict[(penalty1, penalty2)], all_edges_dict, test_set)
        top_20_reults_dict[(penalty1, penalty2, penalty3)] = corr

    return top_20_reults_dict


def run_grid_search2(jacard_index_by_file_dict, all_edges_dict):
    # outlier_panalties = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]
    # bp_qnode_penalties = [0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4]
    outlier_panalties = np.arange(0, 4.2, 0.3)
    bp_qnode_penalties = np.arange(0, 4.2, 0.3)
    qnode_flip_penalties = np.arange(0, 4.2, 2)

    all_penalties_MEM4_dict = {}
    for penalty1 in outlier_panalties:
        for penalty2 in bp_qnode_penalties:
            for penalty3 in qnode_flip_penalties:
                MEM4_by_file_dict = get_all_MEM4_dicts(penalty1, penalty2, penalty3)
                MEM4_by_file_dict = change_keys(MEM4_by_file_dict)
                all_penalties_MEM4_dict[(penalty1, penalty2, penalty3)] = MEM4_by_file_dict

    top_results = get_top_20_results(jacard_index_by_file_dict, all_penalties_MEM4_dict, all_edges_dict)
    num_of_iterations = 5
    for i in range(num_of_iterations-1):
        top_results_iteration = get_top_20_results(jacard_index_by_file_dict, all_penalties_MEM4_dict, all_edges_dict)
        shared_top_results = set(top_results.keys()).intersection(set(top_results_iteration.keys()))
        new_top_results = {}
        for shared_result in shared_top_results:
            new_top_results[shared_result] = top_results[shared_result] + top_results_iteration[shared_result]
        top_results = new_top_results

    for result in top_results.keys():
        top_results[result] = top_results[result]/num_of_iterations

    top_results_sorted = sorted(top_results.items(), key=lambda x: x[1])

    print("\n\n")
    print(top_results_sorted)

    top_20 = [result[0] for result in top_results_sorted]

    # check on 10 test examples
    top_20_reults_dict = {}
    for penalty1, penalty2, penalty3 in top_20:
        corr = print_all_files_correlations_from_dicts2(jacard_index_by_file_dict, all_penalties_MEM4_dict[(penalty1, penalty2, penalty3)], jacard_index_by_file_dict.keys())
        top_20_reults_dict[(penalty1, penalty2, penalty3)] = corr

    print("\nresults for all 20 files:")
    print(top_20_reults_dict)
    return top_results_sorted


def get_csb_set(csb):
    cogs_dict = {}
    for cog in csb:
        if cog not in cogs_dict:
            cogs_dict[cog] = 0
        cogs_dict[cog] = cogs_dict[cog] + 1

    csb_set = {(key, cogs_dict[key]) for key in cogs_dict.keys()}
    return frozenset(csb_set)






def keep_families_from_len_3(csbs_dict):
    new_dict = {}
    for key in csbs_dict:
        csbs_list = csbs_dict[key]
        # if 4 <= len(csbs_list) <= 5:
        if len(csbs_list) >= 3:
            new_dict[key] = csbs_list
    return new_dict


def keep_families_above_50_inst(csbs_dict):
    new_dict = {}
    for key in csbs_dict:
        csbs_list = csbs_dict[key]
        num_instances = int(csbs_list[0].split('\t')[3])
        if num_instances > 50:
            new_dict[key] = csbs_dict[key]

    return new_dict


def keep_families_with_dup(csbs_dict):
    new_dict = {}
    for key in csbs_dict:
        num_dup_cogs = 0    # the number of duplicate cogs
        for (cog, num_cog) in key:
            if num_cog > 1:
                num_dup_cogs += 1

        if num_dup_cogs == 0:
            csbs_list = csbs_dict[key]
            new_dict[key] = csbs_list

    return new_dict




def create_files_from_csbs_dict(csbs_dict):
    first_line = "CSB_ID	Length	Score	Instance_Count	CSB	Main_Category	Family_ID\n"
    families_list = [fam for fam in csbs_dict.values()]
    i=1
    for family in families_list:
        fam_string = first_line[:]
        family_id = str(i) + "-" + (family[0].split('\t'))[6][:-1]
        for csb_line in family:
            fam_string = fam_string + csb_line
        print("\n\n")
        print(family_id)
        print(fam_string)
        i += 1

        path = "all_families_above_4cogs_above50inst_nodup_with_tandem2/" + family_id + ".txt"
        f = open(path, "a")
        f.write(fam_string)
        f.close()


def get_cog_set_abc():
    path = "all_examples/cog_info.txt"
    cog_lines = open(path).readlines()
    cog_set = set()
    for cog_info in cog_lines:
        if "ABC" in cog_info:
            cog_set.add(cog_info.split(";")[0])

    return cog_set


def keep_families_with_abc_cog(csbs_dict):
    cog_set_abc = get_cog_set_abc()
    new_dict = {}
    for key in csbs_dict:
        num_abc_cogs = 0    # the number of ABC-type cogs
        for (cog, num_cog) in key:
            if cog in cog_set_abc:
                num_abc_cogs += 1

        if num_abc_cogs > 0:
            csbs_list = csbs_dict[key]
            new_dict[key] = csbs_list

    return new_dict


def check_tandem(csb):
    new_csb = []
    num_of_tandem_repeats = 0
    for i in range(len(csb)-1):
        if csb[i] != csb[i+1]:
            new_csb.append(csb[i])
        else:
            num_of_tandem_repeats += 1
    new_csb.append(csb[len(csb)-1])

    if num_of_tandem_repeats != 1:
        return csb
    return new_csb


def check_tandem2(csb):
    new_csb = []
    num_of_tandem_repeats = 0
    for i in range(len(csb)-1):
        if csb[i] != csb[i+1]:
            new_csb.append(csb[i])
        else:
            num_of_tandem_repeats += 1
    new_csb.append(csb[len(csb)-1])

    if num_of_tandem_repeats == 0:
        return csb
    return new_csb



def create_new_families_with_duplicatios():
    path = "all_examples/all_families.txt"
    csbs_file_lines = open(path).readlines()

    first_csb = csbs_file_lines[1].split('\t')[4].split(',')
    sum = 0
    csbs_dict = {}
    # csbs_dict_tandem = {}
    for line in csbs_file_lines[1:]:
        line_spitted = line.split('\t')
        csb_len = int(line_spitted[1])
        if 4 <= csb_len <= 9:
        # if 4 <= csb_len <= 5:
        # if csb_len == 3:
            csb = line.split('\t')[4].split(',')
            csb = check_tandem2(csb)
            csb = [cog[:-1] for cog in csb]
            csb_set = get_csb_set(csb)


            if csb_set not in csbs_dict:
                csbs_dict[csb_set] = []

            csbs_dict[csb_set].append(line)
            sum += 1

    print("num of CSBs before", len(csbs_file_lines)-1)

    print("sum", sum)

    print("num of families", len(csbs_dict))

    csbs_dict = keep_families_from_len_3(csbs_dict)
    print("keep_families_from_len_3", len(csbs_dict))

    # csbs_dict = keep_families_with_abc_cog(csbs_dict)
    # print("keep_families_with_abc_cog", len(csbs_dict))

    csbs_dict = keep_families_above_50_inst(csbs_dict)
    print("keep_families_above_50_inst", len(csbs_dict))

    csbs_dict = keep_families_with_dup(csbs_dict)
    print("keep_families_with_dup", len(csbs_dict))


    # create_files_from_csbs_dict(csbs_dict)


def get_csbs_inst_dict():
    path = "all_examples/all_families_instances.fasta"
    file_content = open(path).read()
    content_splitted = file_content.split(">")
    csbs_inst_dict = {}
    for instances in content_splitted:
        csb_id = instances.split("\t")[0]
        # print("\n\n")
        # print(csb_id)
        # print(">" + instances)
        csbs_inst_dict[csb_id] = ">" + instances

    return csbs_inst_dict



def get_fam_inst_string(path, csbs_inst_dict):
    fam_inst_string = ""
    lines = open(path).readlines()

    for line in lines[1:]:
        csb_id = line.split("\t")[0]
        fam_inst_string = fam_inst_string + csbs_inst_dict[csb_id]

    return fam_inst_string


def create_families_instances():
    csbs_inst_dict = get_csbs_inst_dict()
    files_paths = glob.glob(r'all_families_above_4cogs_above50inst_nodup_with_tandem2\*.txt')
    for path in files_paths:
        fam_id = path.split(".")[0]
        fam_path = fam_id + "_instances.fasta"
        print(fam_id)
        print(fam_path)

        fam_inst_string = get_fam_inst_string(path, csbs_inst_dict)
        f = open(fam_path, "a")
        f.write(fam_inst_string)
        f.close()



def keep_families_nested(jacard_index_by_file_dict):
    for file in jacard_index_by_file_dict.keys():
        keys = [key for key in jacard_index_by_file_dict[file].keys()]
        for key in keys[1:]:
            corr_tuple = jacard_index_by_file_dict[file][key]
            if corr_tuple[0] == 0 or corr_tuple[1] == 0 or corr_tuple[2] == 0:
                os.remove("new_families/" + file + ".txt")
                os.remove("new_families/" + file + "_instances.fasta")



def change_families_with_outliers():
    files_paths = glob.glob(r'new_families/*.txt')
    for path in files_paths:
        outlier_index = find_outliers2(path)
        print(outlier_index)

def main():
    start = time.time()

    jacard_index_by_file_dict, all_edges_dict = get_all_jaccard_index_dicts()
    print_all_jaccard_index_dicts(jacard_index_by_file_dict)
    # print(all_edges_dict)


    MEM4_general_by_file_dict, pqtrees_dict = get_all_MEM4_general_dicts(0,1.5, 0, 0, 0, 0, 0, 1)
    MEM4_general_by_file_dict = change_keys(MEM4_general_by_file_dict)
    #
    #
    #5
    breakpoint_by_file_dict = get_all_breakpoint_dicts()
    # # #
    signed_breakpoint_by_file_dict = get_all_signed_breakpoint_dicts()
    #
    #
    # # #
    #
    #
    # #   run correlations and trees diff
    print_all_files_correlations_from_dicts(jacard_index_by_file_dict, MEM4_general_by_file_dict, signed_breakpoint_by_file_dict, pqtrees_dict)

    # print_all_files_correlations_from_dicts_with_edges(jacard_index_by_file_dict, MEM4_general_by_file_dict, all_edges_dict)


    # # trees_diff, trees_diff_AVG= get_trees_difference(jacard_index_by_file_dict, breakpoint_by_file_dict, all_edges_dict)
    # print("\n\ntrees diff: " + str(trees_diff) + "    AVG: " + str(trees_diff_AVG))
    #
    #
    # print_all_files_correlations_from_dicts_by_edges(jacard_index_by_file_dict, MEM4_by_file_dict, all_edges_dict)

    # run_grid_search(jacard_index_by_file_dict, all_edges_dict)


    # top_results = run_grid_search2(jacard_index_by_file_dict, all_edges_dict)


    # run_grid_search3(jacard_index_by_file_dict, signed_breakpoint_by_file_dict)


    # create_new_families_with_duplicatios()
    # # # # # # # #
    # create_families_instances()


    # keep_families_nested(jacard_index_by_file_dict)

    # change_families_with_outliers()

    stop = time.time()
    print("\n\nThe time of the run:", stop - start, "seconds  ;  ", (stop - start)/60, "minutes")

if __name__ == '__main__':
    main()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
