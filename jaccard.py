# Jaccard index calculations

import glob


def calc_nestedness(intersection, instances1, instances2):
    a = intersection
    b = instances1 - intersection
    c = instances2 - intersection
    mini, maxi = min(b, c), max(b, c)
    left = (maxi-mini) / (2*a + mini + maxi)
    # right = a / (a + mini)
    right = 1
    return left * right


def get_jaccard_index(intersection_dict, num_of_instances):
    jaccard_dict = {}
    keys = [k for k in intersection_dict.keys()]
    key1 = keys[0]
    for key2 in intersection_dict:
        jaccard = 1 - (intersection_dict[key1][key2] / (num_of_instances[key1] + num_of_instances[key2] - intersection_dict[key1][key2]))
        sorensen = (num_of_instances[key1] + num_of_instances[key2] - 2*intersection_dict[key1][key2]) / (num_of_instances[key1] + num_of_instances[key2])
        nestedness = calc_nestedness(intersection_dict[key1][key2], num_of_instances[key1], num_of_instances[key2])
        jaccard_dict[key2] = (jaccard, sorensen, nestedness)
    return jaccard_dict


def get_taxa_ranks_dict():
    taxa_ranks = ["kingdom", "phylum", "class", "genus", "species"]
    dict = {}
    for rank in taxa_ranks:
        dict[rank] = {}
    f_taxa_csbfinder = open("taxa_csbfinder.txt")
    for line in f_taxa_csbfinder:
        line_splited = line.split(',')
        uid = (line_splited[0]).split("uid")[1]
        dict["kingdom"][uid] = line_splited[1]
        dict["phylum"][uid] = line_splited[2]
        dict["class"][uid] = line_splited[3]
        dict["genus"][uid] = line_splited[4]
        dict["species"][uid] = line_splited[5][:-1]

    return dict


def get_ranks_dict_of_csbs(path):
    taxa_ranks_dict = get_taxa_ranks_dict()
    taxa_ranks = ["kingdom", "phylum", "class", "genus", "species"]
    f_dataset = open(path)
    ranks_dict_of_csb = {}
    csbID=""
    for line in f_dataset:
        if line[0]=='>':
            csbID = line.split('\t')[0][1:]
            ranks_dict_of_csb[csbID] = {}
            for rank in taxa_ranks:
                ranks_dict_of_csb[csbID][rank] = set()
        else:
            uid = (line.split("uid")[1]).split('\t')[0]
            for rank in taxa_ranks:
                if taxa_ranks_dict[rank][uid] != "-":
                    (ranks_dict_of_csb[csbID][rank]).add(taxa_ranks_dict[rank][uid])

    return ranks_dict_of_csb


def get_num_of_instances_dict(ranks_dict_of_csb):
    taxa_ranks = ["kingdom", "phylum", "class", "genus", "species"]
    num_of_instances_dict = {}

    for rank in taxa_ranks:
        num_of_instances_dict[rank] = {}
        for key in ranks_dict_of_csb:
            num_of_instances_dict[rank][key] = len(ranks_dict_of_csb[key][rank])

    return num_of_instances_dict


def print_jaccard_index_dict(jaccard_index_dict):
    for key1 in jaccard_index_dict.keys():
        line = ""
        for key2 in jaccard_index_dict.keys():
            line = line + str(jaccard_index_dict[key1][key2]) + "  "


def get_intersection_dict(ranks_dict_of_csb):
  taxa_ranks = ["kingdom", "phylum", "class", "genus", "species"]
  intersection_dict = {}

  for rank in taxa_ranks:
      intersection_dict[rank] = {}
      for key1 in ranks_dict_of_csb:
          intersection_dict[rank][key1] = {}
          for key2 in ranks_dict_of_csb:
              ranks_set1 = ranks_dict_of_csb[key1][rank]
              ranks_set2 = ranks_dict_of_csb[key2][rank]
              intersection_set = ranks_set1.intersection(ranks_set2)
              intersection_dict[rank][key1][key2] = len(intersection_set)

  return intersection_dict


def get_jaccard_index_dict(path):
    ranks_dict_of_csb = get_ranks_dict_of_csbs(path)
    intersection_dict = get_intersection_dict(ranks_dict_of_csb)
    num_of_instances = get_num_of_instances_dict(ranks_dict_of_csb)
    lowest_rank = "class"
    jaccard_index_dict = get_jaccard_index(intersection_dict[lowest_rank], num_of_instances[lowest_rank])
    return jaccard_index_dict


def get_all_jaccard_index_dicts():
    print("\n---------- calculating jaccard index ---------------\n\n")
    instances_files_paths = glob.glob(r'input_families/*.fasta')

    jacard_index_by_file_dict = {}
    for path in instances_files_paths:
        file_name = path[15:].split(".")[0]
        jaccard_index_dict = get_jaccard_index_dict(path)
        jacard_index_by_file_dict[file_name] = jaccard_index_dict

    print_all_jaccard_index_dicts(jacard_index_by_file_dict)
    return jacard_index_by_file_dict


def print_all_jaccard_index_dicts(jacard_index_by_file_dict):
    for file_name in jacard_index_by_file_dict:
        print("\n---------------------------------------", file_name, "---------------------------------------\n")
        print("Inverse Jaccard index scores:")
        first_csb = [k for k in jacard_index_by_file_dict[file_name].keys()][0]
        print_dict = {}
        for csb in jacard_index_by_file_dict[file_name].keys():
            print_dict[csb] = jacard_index_by_file_dict[file_name][csb][0]
        print(first_csb, ":     ", print_dict)
