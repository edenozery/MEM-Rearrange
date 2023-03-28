from MEM_general import get_all_MEM_Rearrange_dicts
from signed_breakpoint import get_all_signed_breakpoint_dicts
from MEM import get_all_MEM4_dicts
from jaccard import get_all_jaccard_index_dicts
import time


def main():
    start = time.time()

    bp_qnode_penal, qnode_flip_penal, jump_penal = 1.5, 0.5, 1
    d_T, d_S, delete_T_penal, delete_S_penal = 0, 0, 0, 0

    MEM4_general_by_file_dict = get_all_MEM_Rearrange_dicts(bp_qnode_penal, qnode_flip_penal, d_T, d_S, delete_T_penal, delete_S_penal, jump_penal)

    # MEM_dict = get_all_MEM4_dicts(bp_qnode_penal, qnode_flip_penal, jump_penal)

    # breakpoint_dict = get_all_signed_breakpoint_dicts()

    # jacard_index_by_file_dict = get_all_jaccard_index_dicts()

    stop = time.time()
    print("\n\nThe time of the run:", stop - start, "seconds  ;  ", (stop - start)/60, "minutes")


if __name__ == '__main__':
    main()


