import numpy as np
import time
import sys
start = time.clock()  # record running time


def reg2hicrep(chr, res, start_bin, end_bin, reg_out, hicrep_in):
    # compute dimension of HiC matrix
    dim = round((int(end_bin) - int(start_bin))/int(res) + 1)

    # initialize zero matrix
    initial_matrix = np.zeros((dim, dim), dtype=int)
    matrix_listed = initial_matrix.tolist()

    # generate bins
    bins = []
    i = 0
    while i < dim:
        bins.append(np.arange(i * int(res), i * int(res) + int(res) + 1, int(res)).tolist())
        i += 1

    with open(reg_out) as infile:
        infile_read = infile.readlines()[1:]
        # fill zero matrix with contact number
        for line in infile_read:
            sline = line.strip().split()
            if sline[0] == "bin1_mid":
                continue
            bin1 = round((int(float(sline[1])) - int(start_bin))/int(res))
            bin2 = round((int(float(sline[2])) - int(start_bin))/int(res))
            interaction = int(sline[3])
            matrix_listed[bin1][bin2] = interaction
        with open(hicrep_in, "w+") as outfile:
            for bin, num in zip(bins, matrix_listed):
                outfile.write("chr" + str(chr) + "\t" + "\t".join(str(j) for j in bin) + "\t" + "\t".join(str(i) for i in num) + "\n")

    end = time.clock()
    print(end - start)


if __name__ == "__main__":
    if len(sys.argv) != 7:
        print("usage : reg_to_hicrep")
        sys.exit()

    reg2hicrep(chr=sys.argv[1], res=sys.argv[2], start_bin=sys.argv[3], end_bin=sys.argv[4],
               reg_out=sys.argv[5], hicrep_in=sys.argv[6])
