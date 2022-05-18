# Break point distance calculations

# calculating break-point distance

import glob

def get_bp_from_csbs(csb1, csb2):
    csb1_list = [csb[:-1] for csb in csb1.split(',')]
    csb2_list = [csb[:-1] for csb in csb2.split(',')]
    bp = 0
    for i in range(1, len(csb1_list)):
        set1 = set([csb1_list[i - 1], csb1_list[i]])
        found = False
        for j in range(1, len(csb2_list)):
            set2 = set([csb2_list[j - 1], csb2_list[j]])
            if (set1 == set2):
                found = True

        if (not found):
            bp += 1
    return bp

def change_sign(sign):
    if sign == "+":
        return "-"
    return "+"

def get_bp_from_children_orders(csb1, csb2):
  co_tree = [(csb[:-1], csb[-1]) for csb in csb1.split(',')]
  co_string = [(csb[:-1], csb[-1]) for csb in csb2.split(',')]
  # print("co_tree", co_tree)
  # print("co_string", co_string)
  if len(co_string) > len(co_tree):
      temp = co_string
      co_string = co_tree
      co_tree = temp
  bp = 0
  for i in range(1,len(co_tree)):
    tup1 = (co_tree[i - 1][0], co_tree[i][0])
    found = False
    for j in range(1,len(co_string)):
      tup2 = (co_string[j - 1][0], co_string[j][0])
      if tup1 == tup2 and co_tree[i - 1][1] == co_string[j - 1][1] and co_tree[i][1] == co_string[j][1]:
        found = True
        break

      tup2 = (co_string[j][0], co_string[j-1][0])
      if tup1 == tup2 and co_tree[i][1] == change_sign(co_string[j - 1][1]) and co_tree[i-1][1] == change_sign(co_string[j][1]):
        found = True
        break


    if(not found):
      bp += 1

  # #print(bp)
  return bp



def get_breakpoint_from_file(csbs_file_path):
    csbs_file_lines = open(csbs_file_path).readlines()
    csbs_file2_lines = open(csbs_file_path).readlines()
    bp_dict = {}
    for line in csbs_file_lines:
        if line[0] != 'C' and line[0] != "":
            splited_line = line.split('\t')
            bp_dict[splited_line[0]] = {}
            csb1 = splited_line[4]
            # print(csb1)

            i = 0
            for line2 in csbs_file2_lines:
                i += 1
                if line2[0] != 'C':
                    splited_line2 = line2.split('\t')
                    csb2 = splited_line2[4]
                    bp_dict[splited_line[0]][splited_line2[0]] = get_bp_from_children_orders(csb1, csb2)

    return bp_dict


def print_dict(dict):
    for csb in dict.keys():
        print(dict[csb])

def print_all_jaccard_index_dicts(jacard_index_by_file_dict):
    for key in jacard_index_by_file_dict:
        print("\n\n----------- ", key, " -------------")
        print_dict(jacard_index_by_file_dict[key])


def get_all_signed_breakpoint_dicts():
    print("\n\n---------- calculating break point ---------------\n")
    csbs_files_paths = glob.glob(r'new_families\*.txt')

    breakpoint_by_file_dict = {}
    for path in csbs_files_paths:
        print(path)
        breakpoint_dict = get_breakpoint_from_file(path)
        file_name = path[13:].split(".")[0]
        breakpoint_by_file_dict[file_name] = breakpoint_dict

    return breakpoint_by_file_dict


def main():
    # breakpoint_by_file_dict = get_all_breakpoint_dicts()
    # print_all_jaccard_index_dicts(breakpoint_by_file_dict)
    print(get_breakpoint_from_file(r'all_examples/14796.txt'))

    # print("------------------------------------------------------------")
    # print(get_bp_from_children_orders("COG1117-,COG0581-,COG0573-,COG0226-,COG0704-", "COG1117-,COG0581-,COG0573-,COG0226-,COG0704-"))


if __name__ == '__main__':
    main()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
