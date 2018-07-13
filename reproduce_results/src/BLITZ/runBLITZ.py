from SPPformat_to_csc import SPPformat_to_csc
import scipy.sparse as sp
import numpy as np
import time
import blitzl1
import csv
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--results_dir", type=str,
                    help="Directory where to write the results.",
                    default="../results/Blitz")

parser.add_argument("--data_file", type=str,
                    help="The file that contains the data, in SPP format.",
                    default="../data_preprocessed/Bernoulli_n1000_p1000_" +
                            "qunif_coefnormal_rs0_nnzd.tsv")

parser.add_argument("--nlambda", type=int,
                    help="The number of values of the regularisation" +
                    "parameter lambda for which to compute a solution",
                    default=100)

parser.add_argument("--lambdaMinRatio", type=float,
                    help="The ratio between the largest lambda (the one for" +
                    "which one feature enters the regularisation path), and" +
                    "the smallest value of lambda tested", default=0.01)

parser.add_argument("--maxSelectedFeatures", type=int,
                    help="The maximum number of features that can be" +
                    "selected before the algorithm is stopped", default=150)

parser.add_argument("--useBias", type=int,
                    help="whether or not to use a bias", default=1)

parser.add_argument("--tol", type=float,
                    help="The tolerance for which Blitz should stop",
                    default=1)

args = parser.parse_args()

# Save results in files
pos_start = args.data_file.rfind("/")
pos_end = args.data_file.rfind(".")
filename = args.data_file[pos_start+1:pos_end]
filename_spec = (filename + "_nlam" + str(args.nlambda) + "_rat" +
                 str(args.lambdaMinRatio) +
                 "_max" + str(args.maxSelectedFeatures) + "_Bi" +
                 str(args.useBias) + "_tol" + str(args.tol))

file_model = open(args.results_dir + "/" + filename_spec + "_model.csv", "wb")
csv_writer_model = csv.writer(file_model, delimiter=",")

file_preproc = open(args.results_dir + "/" + filename_spec +
                    "_preproc.csv", "wb")
csv_writer_preproc = csv.writer(file_preproc, delimiter=",")
csv_writer_preproc.writerow("time_preproc")

file_stats = open(args.results_dir + "/" + filename_spec + "_stats.csv", "wb")
csv_writer_stats = csv.writer(file_stats, delimiter=",")
csv_writer_stats.writerow(["lambda", "n_features_selected", "time"])

# Load the binary mutation matrix
data = SPPformat_to_csc(file=args.data_file)
n, p = data['X'].shape

# Center the response
y = data['y'] - sum(data['y'])/n

# Create the sparse matrix that will hold interactions
indices = []
indptr = np.empty(shape=p*(p+1)/2+1, dtype=int)
names = np.empty(shape=p*(p+1)/2, dtype=tuple)
tot = 0
start = 0
end = p
for j in range(p):
    diag = sp.diags(data['X'][:, j].toarray().reshape((n,)),
                    format='csc', dtype=int)
    Xint_sub = diag.dot(data['X'][:, j:p])
    indices.append(Xint_sub.indices)
    tmp_indptr = Xint_sub.indptr
    indptr[start:(end+1)] = tmp_indptr + tot
    names[start:end] = [(j, k) for k in range(j, p)]
    tot += Xint_sub.getnnz()
    start = end
    end += Xint_sub.shape[1] - 1

indices = np.hstack(indices)
indices = indices.astype('int64')
Xint = sp.csc_matrix((np.ones(tot, dtype=np.int8),
                     indices, indptr), dtype=np.int8, shape=(n, p*(p+1)/2))


# Run blitz on SNPs data with interactions
blitzl1.set_tolerance(args.tol)
blitzl1.set_verbose(False)
blitzl1.set_use_intercept(args.useBias)
prob = blitzl1.LassoProblem(Xint, y)

# Compute lambda_max
t0_lammax = time.time()
lammax = prob.compute_lambda_max()
t1_lammax = time.time()
csv_writer_preproc.writerow([t1_lammax - t0_lammax])
file_preproc.close()

# Define the values of lambda for which the solution will be computed
lam = [lammax*pow(10, np.log10(args.lambdaMinRatio)*t/args.nlambda)
       for t in range(1, args.nlambda+1)]

# Solve the LASSO for every value of lambda in the regularisation path
initial_x = np.zeros(Xint.shape[1])
initial_intercept = 0
for i in range(0, args.nlambda):
    print i
    t0 = time.time()
    sol = prob.solve(lam[i],
                     initial_x=initial_x,
                     initial_intercept=initial_intercept)
    initial_x = sol.x
    initial_intercept = sol.intercept
    t1 = time.time()

    ind_nnz = np.nonzero(sol.x)[0]
    csv_writer_stats.writerow([lam[i], len(ind_nnz), t1-t0])
    for _, k in enumerate(ind_nnz):
        csv_writer_model.writerow([i, lam[i], names[k][0], names[k][1],
                                   sol.x[k]])

    # print i, lam[i], len(np.nonzero(sol.x)[0])
    if len(ind_nnz) >= args.maxSelectedFeatures:
        break

file_model.close()
file_stats.close()
