# WHInter
WHInter allows you to efficiently fit a LASSO to binary matrices with up to a million features, which corresponds to O(10<sup>12</sup>) interactions. For example, on a problem with 1000 samples and 10.000 features (corresponding to 10<sup>10</sup> interactions), WHINter fits a LASSO for a whole regularisation path in less than one minute. WHInter is described in the paper [WHInter: A Working set algorithm for High-dimensional sparse second order Interaction models] (http://proceedings.mlr.press/v80/morvan18a/morvan18a.pdf).

#### How to run WHINter
##### Environment
WHINter is entirely written in C++. It was tested with:
* gcc version 4.9.4
* GNU Make 3.81

##### Compilation
To compile:
```bash
make
```

##### Command line arguments
To use WHInter:
```bash
./src/train_WHInter [options] [filename]
```
where options are:
* `-nlambda`: The number of values of the regularisation parameter lambda for which to compute a solution (default 100).
* `-lambdaMinRatio`: The ratio between the largest lambda (the one for which one feature enters the regularisation path), and the smallest value of lambda tested (default 0.01).
* `-maxSelectedFeatures`: The maximum number of features that can be selected before the algorithm is stopped (default 150).
* `-useBias`: Use a bias (1) ot not (0) (default 1).
* `-useMyMips`: The type of MIPS solver to be used:
    * 0: Naive solver
    * 1: Experimental solver
    * 2: IL (inverted index)
The default is 2. WHInter is presented in the paper with option 2, which is by far the fastest option.
* `-typeBound`: The type of branch bound eta to be used:
    * 0: alpha is taken equal to 1 (eta_1).
    * 1: alpha is the minimiser of eta (eta_{min}).
    * 2: (eta_{alpha_{\ell_2}}).
The default is 2.
* `-F`: The subproblem solver performs dynamic screening every F iterations (default 50).
* `-eps`: Convergence threshold on the relative duality gap for the subproblem solver (default 1e-8).
* `-pathResults`: Directory where to save the results (default "./")

Example:
```bash
mkdir -p ./results/
./src/train_WHInter -pathResults ./results/ ./toy_data/data_n1000_p1000.tsv
```

#### Input file fromat
WHInter takes as input text files such as [this toy example](./toy_data/data_n1000_p1000.tsv). Each line corresponds to a sample. The first element of the line is the response value, and the following elements describe the feature values in the [index]:[value] sparse format, where index is the feature index and value is the corresponding non-zero value.

If needed, the Python function `binary_to_sparse_format` defined in [simulate_lasso.py](./reproduce_results/get_data_preprocessed/simulate_lasso.py) allows to convert a design matrix and a response vector available as numpy arrays or pandas dataframes format into the right format.

#### Output files
There are 3 output files for each run ending in `_preproc.csv`, `_stat.csv` and `_model.csv` respectively. These files describe the learned model for every regularisation parameter tested as well as information about the run.

###### `*_model.csv`
Records the feature IDs that were selected for each regularisation parameter as well as their corresponding weight in the model. Feature IDs are reported as pairs (`Feature_ID1`, `Feature_ID2`) where (5,10) is for example the interaction feature between the 5<sup>th</sup> and 10<sup>th</sup> columns of the original matrix. In particular, the feature (5, 5) simply represents the 5<sup>th</sup> feature in the orginal matrix. The columns in this file represent, in that order:
* the index of the regularisation paramter
* the value of the regularisation parameter
* `Feature_ID1`
* `Feature_ID2`
* corresponding coefficient in the optimal solution.

The set of selected features corresponding to a given regularisation parameter can therefore be found by taking all the lines with the same regularisation parameter.

###### `*_preproc.csv`
Records the time (in milliseconds) to find \lambda_{max} and to initialise m^{ref} (line 5 in Algorithm 2).

###### `*_stats.csv`
Each line of this file corresponds to one regularisation parameter tested. Each column contains the following information:
* `lambda`: The current regularisation parameter.
* `n_features_selected`: Number of features with non-zero weights at the optimum for a given regularisation parameter.
* `time`: Total time to solve the problem for a given regularisation parameter.
* `iterB`: Number of outer loops (line 10 in Algorithm 2) performed for a given regularisation parameter.
* `n_fictive_active`: Number of branches which could not have been pruned based on the branch bound zeta for the current dual point.

The remaining columns describe quantities that change for each outer loop, for example the number of branches that can be pruned. For this reason we use for each quantity as many columns as there are outer loops. For example the column named `col` becomes `col1` and `col2` if there are 2 outer loops. Note that the number of outer loops performed may differ according to the regularisation parameter used. Therefore the maximum number of outer loops across regularisation parameters is used as a reference and for regularisation parameters for which less outer loops were perfomed we fill the corresponding columns with zeros. 
* `n_branch_opened`: Number of branches that could not be pruned.
* `n_activeset`: Working set size.
* `nb_iter_CD`: Number of iterations of the inner solver.
* `time_check_branch`: Time (in milliseconds) for computing alpha, computing the value of the branch bound, and checking wether the branch can be pruned. This corresponds to lines 11 to 17 in Algorithm 2.
* `time_MIPS`: Time (in milliseconds) to perform maximum inner product search (line 18 in Algorithm 2)
* `time_CD`: Time (in milliseconds) to solve the subproblem restricted to the working set (line 24 in Algorithm 2).

#### Reproducing the results of the paper
In this project, Python (version 2.7.10) and R (version 3.2.4) were used to simulate data, parse WHInter outputs and plot figures, or run other methods we compare with. Python was used with the following modules:
* scipy 0.17.0
* numpy 1.11.0
* csv 1.0
* argparse 1.1
* pandas 0.18.0
* matplotlib 1.5.1
* time
* blitzl1 (see BLITZ [github page](https://github.com/tbjohns/BlitzL1))

To reproduce the simulation results of the [paper](http://proceedings.mlr.press/v80/morvan18a/morvan18a.pdf):
```bash
./reproduce_results/reproduce_simulations.sh
```
This script generates simulated data, computes the LASSO with second-order interaction terms with WHInter and competing methods (this takes a few hours) and plot figures based on the obtained results.

#### Upcoming updates
For now, WHInter only takes binary data as input. We plan to improve the code so that it can take as input sparse data in general. Moreover, only the LASSO is available for the moment. We plan to make sparse logistic regression available for classification tasks.

Any feedback on this early version is greatly appreciated :)




