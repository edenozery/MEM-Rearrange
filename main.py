from MEM_general import *
import time


def main():
    start = time.time()

    MEM4_general_by_file_dict, pqtrees_dict = get_all_MEM_Rearrange_dicts(0,1.5, 0.5, 0, 0, 0, 0, 1)

    stop = time.time()
    print("\n\nThe time of the run:", stop - start, "seconds  ;  ", (stop - start)/60, "minutes")


if __name__ == '__main__':
    main()


