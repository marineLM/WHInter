// This code was taken from https://github.com/takeuchi-lab/SafePatternPruning/tree/master/itemLasso 
// This file has been modified compared to the original one:
//  - Added get_maxval2 and safeprune2 functions
//  - Changed the calls get_maxval and safeprune to get_maxval2 and safeprune2
//  - Modified the function "read" which reads command line arguments 
//  - Modified how the sequence of lambda is created for more flexibility
//  - Modified the variable names corresponding to the command line arguments
//  - Added the condition if (model.size() >= maxSelectedFeatures) break; so that
//  computation of the regularization path is stopped if a certain number of features
//  has entered the model

#include <cstdio>
#include <cstdlib>
#include <random>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>
#include <ctime>
#include <string>
#include <string.h>
#include <armadillo>

using namespace std;

struct node {
	vector<int> key;
	vector<int> x;
	double val;
};

struct feature {
	vector<int> x;
	double w;
};

struct ball {
	vector<double> center;
	double radius;
};

typedef map< vector<int>, feature > model;

int nlambda;
double lambdaMinRatio;
int maxSelectedFeatures;
int useBias;
int F;
double eps;
int maxdepth = 2;
string pathResults;

int n, d;

vector< vector<int> > Z;
vector<double> Y;

void read(int argc, char **argv);
void read_data(char *filename);
node get_maxval(const vector<double> &v);
void get_maxval(const vector<double> &v, vector<int> &parentkey, const vector<int> &parentx, int depth, int j, node &maxnode, int n_branch_opened);
node get_maxval2(const vector<double> &v, int &n_branch_opened);
void safeprune(const ball &ball, vector<node> &safeset);
void safeprune(const ball &ball, vector<node> &safeset, vector<int> &parentkey, const vector<int> &parentx, int depth, int j);
void safeprune2(const ball &ball, vector<node> &safeset, int &n_branch_opened);

int main(int argc, char **argv) {

	read(argc, argv);

	string filename = argv[argc - 1];
    size_t pos_start = filename.find_last_of("/");
    size_t pos_end = filename.find_last_of(".");
    string file = filename.substr(pos_start + 1, pos_end - (pos_start + 1));
	
	string str_rat = to_string(lambdaMinRatio);
    str_rat.erase(str_rat.find_last_not_of('0') + 1, string::npos);
    char str_eps [10];
    sprintf(str_eps, "%g", eps); 
	ofstream myfile;
	myfile.open(pathResults +
                file + "_nlam" + to_string(nlambda) + "_rat" +  str_rat + 
                "_max" + to_string(maxSelectedFeatures) + "_Bi" + to_string(useBias) +
                "_F" + to_string(F) + "_eps" + str_eps + "_stats.csv");
	ofstream file_submodel_size;
	file_submodel_size.open(pathResults +
                file + "_nlam" + to_string(nlambda) + "_rat" +  str_rat + 
                "_max" + to_string(maxSelectedFeatures) + "_Bi" + to_string(useBias) +
                "_F" + to_string(F) + "_eps" + str_eps + "_cd.csv");
	ofstream file_model;
	file_model.open(pathResults +
                file + "_nlam" + to_string(nlambda) + "_rat" +  str_rat + 
                "_max" + to_string(maxSelectedFeatures) + "_Bi" + to_string(useBias) +
                "_F" + to_string(F) + "_eps" + str_eps + "_model.csv");
	ofstream file_preproc;
	file_preproc.open(pathResults +
                file + "_nlam" + to_string(nlambda) + "_rat" +  str_rat + 
                "_max" + to_string(maxSelectedFeatures) + "_Bi" + to_string(useBias) +
                "_F" + to_string(F) + "_eps" + str_eps + "_preproc.csv");

	auto t_start = std::chrono::high_resolution_clock::now();

	double bias = 0;

	// residual
	vector<double> r(n);
	for (int i = 0; i < n; i++) r[i] = Y[i];

	file_preproc << "time_preproc, n_branch_opened" << endl;
	int tmp = 0;
	auto t_preproc_start = std::chrono::high_resolution_clock::now();
	node maxnode = get_maxval2(r, tmp);
	auto t_preproc_end = std::chrono::high_resolution_clock::now();
	file_preproc << std::chrono::duration_cast<std::chrono::milliseconds>(t_preproc_end - t_preproc_start).count() << ", " << tmp << endl;
	double lammax = maxnode.val;
	cout << "------------------------------" << endl;
	cout << "lambda_max = " << lammax << endl;

	model model;
	model[maxnode.key] = feature{maxnode.x, 0};


	vector<double> theta(n);
	for (int i = 0; i < n; i++) theta[i] = r[i]/lammax;

	// vector<double> lambda(T-1);
	// for (int t = 0; t < T-1; t++) {
	// 	lambda[t] = lammax*pow(10,-2.0*(t+1)/(T-1));
	// }

	vector<double> lambda(nlambda);
    for (int t = 1; t <= nlambda; t++) {
        lambda[t-1] = lammax*pow(10,log10(lambdaMinRatio)*t/nlambda);
    }

	//
	// START: solution path algorithm
	//
	myfile << "lambda" << " " << "n_branch_opened_proj" << " ";
	myfile << "n_branch_opened_safepruning" << " " << "n_safeset" << " " << "n_features_selected" << " ";
	myfile << "nb_iter_CD" <<  " " << "time" << " " << "time_proj1" << " " << "time_safeprune" << " " << "time_CD" << "\n";
	for (int t = 0; t < nlambda; t++) {
	// for (int t = 1; t <= 60; t++) {
		auto t0 = std::chrono::high_resolution_clock::now();
		std::chrono::high_resolution_clock::time_point t_start_proj1 = std::chrono::high_resolution_clock::now();
		std::chrono::high_resolution_clock::time_point t_end_proj1 = std::chrono::high_resolution_clock::now();
		std::chrono::high_resolution_clock::time_point t_start_safeprune = std::chrono::high_resolution_clock::now();
		std::chrono::high_resolution_clock::time_point t_end_safeprune = std::chrono::high_resolution_clock::now();
		std::chrono::high_resolution_clock::time_point t_start_CD = std::chrono::high_resolution_clock::now();
		std::chrono::high_resolution_clock::time_point t_end_CD = std::chrono::high_resolution_clock::now();

		double lam = lambda[t];
		cout << "------------------------------" << endl;
		printf("[%d] lambda: %.9f (%.9f)\n", t, lam, log10(lam/lammax));
		myfile << lam << " ";

		// warm start
		int m = model.size();
		vector<int> index(m);
		vector<double> w(m);
		vector< vector<int> > x(m);
		vector< vector<int> > key(m);
		int j = 0;
		for (auto it = model.begin(); it != model.end(); it++) {
			index[j] = j;
			w[j] = it->second.w;
			x[j] = it->second.x;
			key[j] = it->first;
			j++;
		}

		// START: pre_solve
		double P_old = 10000000;
		for (int iter = 0; iter <= 1000000; iter++) {

			for (int j = 0; j < m; j++) {
				int i = j + rand()%(m-j);
				swap(index[i], index[j]);
			}

			double L1norm = 0;
			for (int s = 0; s < m; s++) {
				int j = index[s];

				double xTr = 0;
				for (int k = 0; k < (int)x[j].size(); k++) {
					int idx = x[j][k];
					r[idx] += w[j];
					xTr += r[idx];
				}

				if (xTr > lam) {
					w[j] = (xTr-lam)/x[j].size();
				} else if (xTr < -lam) {
					w[j] = (xTr+lam)/x[j].size();
				} else {
					w[j] = 0;
				}

				for (int k = 0; k < (int)x[j].size(); k++) {
					int idx = x[j][k];
					r[idx] -= w[j];
				}
				L1norm += fabs(w[j]);
			}

			// bias
			if (useBias) {
				double tmp = 0;
				for (int i = 0; i < n; i++) {
					r[i] += bias;
					tmp += r[i];
				}
				bias = tmp/n;
				for (int i = 0; i < n; i++) {
					r[i] -= bias;
				}
			}

			double loss = 0;
			for (int i = 0; i < n; i++) {
				loss += r[i]*r[i];
			}
			double P_new = 0.5*loss + lam*L1norm;
			if (fabs((P_old-P_new)/P_old) < 1e-8) break;
			P_old = P_new;
		}
		// END: pre_solve

		double D = 0;
		for (int i = 0; i < n; i++) {
			D += 0.5*Y[i]*Y[i] - 0.5*lam*lam*(theta[i]-Y[i]/lam)*(theta[i]-Y[i]/lam);
		}

		//
		// START: coodinate descent
		//
		for (int iter = 1; iter <= 1000000; iter++) {
			// START: dual
			if (iter % F == 1) {
				double loss = 0;
				double yTr = 0;
				for (int i = 0; i < n; i++) {
					loss += r[i]*r[i];
					yTr += Y[i]*r[i];
				}

				// L1norm
				double L1norm = 0;
				for (int s = 0; s < m; s++) {
					int j = index[s];
					L1norm += fabs(w[j]);
				}

				double maxval = 0;
				if (iter == 1) {
					t_start_proj1 = std::chrono::high_resolution_clock::now();
					int n_branch_opened_proj = 0;
					node maxnode = get_maxval2(r, n_branch_opened_proj);
					myfile << n_branch_opened_proj << " ";
					if (maxnode.key.size() == 2){
						cout << "maxval found: " << maxnode.key[0] <<  " " << maxnode.key[1] << endl;
					}else{
						cout << "maxval found: " << maxnode.key[0] << endl;
					}
					maxval = maxnode.val;
					t_end_proj1 = std::chrono::high_resolution_clock::now();
				} else {
					for (int s = 0; s < m; s++) {
						int j = index[s];
						int xnorm = x[j].size();
						double xtr = 0;
						for (int i = 0; i < xnorm; i++) {
							int idx = x[j][i];
							xtr += r[idx];
						}
						if (fabs(xtr) > maxval) maxval = fabs(xtr);
					}
				}

				// dual feasible solution
				double alpha = min(max(yTr/(lam*loss), -1/maxval), 1/maxval);

				double P = 0.5*loss + lam*+L1norm;
				double D_new = -0.5*lam*lam*alpha*alpha*loss + lam*alpha*yTr;
				if (D < D_new) {
					D = D_new;
					for (int i = 0; i < n; i++) {
						theta[i] = alpha*r[i];
					}
				}

				double gap = P-D;
				if (gap/P < eps) {
					// update model
					cout << "loss: " << loss << " L1norm: " << L1norm << endl;
					model.clear();
					int active = 0;
					for (int s = 0; s < m; s++) {
						int j = index[s];
						if (w[j] != 0) {
							active++;
							model[key[j]] = feature{x[j], w[j]};
						}
					}
					printf("[iter %4d] primal: %.9f, dual: %.9f, gap(relative): %.9f, active: %d\n", iter, P, D, (P-D)/P, active);
					file_submodel_size << active << endl;
					t_end_CD = std::chrono::high_resolution_clock::now();
					auto t1 = std::chrono::high_resolution_clock::now();
					if (iter == 1) {myfile << 0 << " " << 0 << " ";}
					myfile << active << " ";
					myfile << iter << " ";
					myfile << std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count() << " ";
					myfile << std::chrono::duration_cast<std::chrono::milliseconds>(t_end_proj1 - t_start_proj1).count() << " ";
					myfile << std::chrono::duration_cast<std::chrono::milliseconds>(t_end_safeprune - t_start_safeprune).count() << " ";
					if (iter==1){
						myfile << 0 << endl;
					}else{
						myfile << std::chrono::duration_cast<std::chrono::milliseconds>(t_end_CD - t_start_CD).count() << endl;
					}
					cout << "total time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t1-t0).count() << endl;
					for (auto it = model.begin(); it != model.end(); it++) {
						vector<int> vec = it->first;
						if (vec.size() == 2){
							cout << vec[0] << " " << vec[1] << " :" << it->second.w << endl;
							file_model << t << ", " << lam << ", " << vec[0] << ", " << vec[1] << ", " << it->second.w << endl;
						}else{
							cout << vec[0] << " :" << it->second.w << endl;
							file_model << t << ", " << lam << ", " << vec[0] << ", " << vec[0] << ", " << it->second.w << endl;
						}
					}
					break;
				}

				// START: safe pruning
				double radius = sqrt(2*gap)/lam;
				vector<node> safeset;
				int active = 0;
				if (iter == 1) { // first pruning
					t_start_safeprune = std::chrono::high_resolution_clock::now();
					map<vector<int>, feature> old;
					for (int j = 0; j < m; j++) {
						if (w[j] != 0) {
							old[key[j]] = feature{x[j], w[j]};
						}
					}

					ball B = {theta, radius};
					int n_branch_opened_sp = 0;
					safeprune2(B, safeset, n_branch_opened_sp);
					myfile << n_branch_opened_sp << " ";

					m = (int)safeset.size();
					myfile << m << " ";
					index.resize(m);
					w.resize(m);
					x.resize(m);
					key.resize(m);
					for (int j = 0; j < m; j++) {
						auto flag = old.find(safeset[j].key);
						if (flag != old.end()) {
							w[j] = flag->second.w;
						} else {
							w[j] = 0;
						}
						if (w[j] != 0) active++;
						x[j] = safeset[j].x;
						key[j] = safeset[j].key;
						index[j] = j;
					}
					t_end_safeprune = std::chrono::high_resolution_clock::now();
				} else { // dynamic screening
					for (int s = 0; s < m; s++) {
						int j = index[s];
						int xnorm = x[j].size();
						double xtc = 0;
						for (int i = 0; i < xnorm; i++) {
							int idx = x[j][i];
							xtc += theta[idx];
						}
						if (w[j] != 0) active++;
						if (fabs(xtc)+radius*sqrt(xnorm) < 1) {
							w[j] = 0;
							m--;
							swap(index[s], index[m]);
							s--;
						}
					}
				}
				printf("[iter %4d] primal: %.9f, dual: %.9f, gap(relative): %.9f, safeset: %d, active: %d\n", iter, P, D, (P-D)/P, m, active);
				file_submodel_size << m << " ";
			}
			// END: dual

			if (iter == 1) t_start_CD = std::chrono::high_resolution_clock::now();
			for (int j = 0; j < m; j++) {
				int i = j + rand()%(m-j);
				swap(index[i], index[j]);
			}

			// START: update wj
			for (int s = 0; s < m; s++) {
				int j = index[s];

				double xTr = 0;
				for (int k = 0; k < (int)x[j].size(); k++) {
					int idx = x[j][k];
					r[idx] += w[j];
					xTr += r[idx];
				}

				if (xTr > lam) {
					w[j] = (xTr-lam)/x[j].size();
				} else if (xTr < -lam) {
					w[j] = (xTr+lam)/x[j].size();
				} else {
					w[j] = 0;
				}

				for (int k = 0; k < (int)x[j].size(); k++) {
					int idx = x[j][k];
					r[idx] -= w[j];
				}
			}
			// END: update wj

			// bias
			if (useBias) {
				double tmp = 0;
				for (int i = 0; i < n; i++) {
					r[i] += bias;
					tmp += r[i];
				}
				bias = tmp/n;
				for (int i = 0; i < n; i++) {
					r[i] -= bias;
				}
			}
		}

		if (model.size() >= maxSelectedFeatures) break;
		//
		// END: coodinate descent
		//

	}
	//
	// END: solution path algorithm
	//
	auto t_end = std::chrono::high_resolution_clock::now();
	cout << "total runtime: " << std::chrono::duration_cast<std::chrono::milliseconds>(t_end-t_start).count() << endl;
	myfile.close();
	file_submodel_size.close();
	file_model.close();
	file_preproc.close();
	return 0;
}

node get_maxval(const vector<double> &v) {
	node maxnode;
	maxnode.val = 0;

	vector<int> key;
	key.reserve(maxdepth);
	get_maxval(v, key, vector<int>(), 0, 0, maxnode, 0);
	return maxnode;
}

void get_maxval(const vector<double> &v, vector<int> &parentkey, const vector<int> &parentx, int depth, int j, node &maxnode, int n_branch_opened) {
	for (; j < d; j++) {
		parentkey.push_back(j);

		vector<int> x; // and
		if (depth) {
			set_intersection(Z[j].begin(), Z[j].end(), parentx.begin(), parentx.end(), back_inserter(x));
		} else {
			x = Z[j];
		}

		if (!x.size()) {
			parentkey.pop_back();
			continue;
		}

		// ここから内積計算
		double xTv = 0;
		double p = 0;
		double m = 0;
		for (int i = 0; i < (int)x.size(); i++) {
			double val = v[x[i]];
			xTv += val;
			(val > 0) ? p += val : m += val;
		}

		if (max(p,-m) < maxnode.val) {
			parentkey.pop_back();
			continue;
		}

		double score = fabs(xTv);
		if (score > maxnode.val) {
			maxnode.key = parentkey;
			maxnode.x = x;
			maxnode.val = score;
		}

		if (j < d-1 && depth+1 < maxdepth) {
			if (depth == 0){ n_branch_opened ++; }
			get_maxval(v, parentkey, x, depth+1, j+1, maxnode, n_branch_opened);
		}

		parentkey.pop_back();
	}
	if (depth == 0){ cout << "# branch opened: " << n_branch_opened << endl;}
	// if (depth == 0){ myfile << n_branch_opened << " ";}
}

node get_maxval2(const vector<double> &v, int &n_branch_opened){
    node maxnode;
    maxnode.val = 0;
    double score;
    vector<int> branch_not_screened;
    // Here
    branch_not_screened.reserve(d);
    for (int i = 0; i < d; i++){
        vector<int> x = Z[i];
        double xTv = 0;
        double p = 0;
        double m = 0;
        for (int k = 0; k < (int)x.size(); k++) {
            double val = v[x[k]];
            xTv += val;
            (val > 0) ? p += val : m += val;
        }
        if (max(p,-m) < maxnode.val) {
            continue;
        }
        branch_not_screened.push_back(i);
        score = fabs(xTv);
        if (score > maxnode.val) {
            maxnode.key = vector<int>{i};
            maxnode.x = x;
            maxnode.val = score;
        }
        for (int b = 0; b < branch_not_screened.size(); b++){
            int j = branch_not_screened[b];
            if (i==j) continue;
            vector<int> xx; //and
            set_intersection(x.begin(), x.end(), Z[j].begin(), Z[j].end(), back_inserter(xx));
            double xxTv = 0;
            for (int k = 0; k < xx.size(); k++){
                xxTv += v[xx[k]];
            }
            score = fabs(xxTv);
            if (score > maxnode.val) {
                maxnode.key = vector<int>{i, j};
                maxnode.x = xx;
                maxnode.val = score;
            }
        }
    }
	n_branch_opened = branch_not_screened.size();
	cout << "# branch opened: " << n_branch_opened << endl;
    return maxnode;
}

void safeprune(const ball &ball, vector<node> &safeset) {
	vector<int> key;
	key.reserve(maxdepth);
	safeprune(ball, safeset, key, vector<int>(), 0, 0);
}

void safeprune(const ball &ball, vector<node> &safeset, vector<int> &parentkey, const vector<int> &parentx, int depth, int j) {
	for (; j < d; j++) {
		parentkey.push_back(j);

		vector<int> x; //and
		if (depth) {
			set_intersection(Z[j].begin(), Z[j].end(), parentx.begin(), parentx.end(), back_inserter(x));
		} else {
			x = Z[j];
		}

		if (!x.size()) {
			parentkey.pop_back();
			continue;
		}

		// centerとxの内積
		double xTc = 0;
		double p = 0;
		double m = 0;
		for (int i = 0; i < (int)x.size(); i++) {
			double val = ball.center[x[i]];
			xTc += val;
			(val > 0) ? p += val : m += val;
		}

		if (max(p, -m) + ball.radius*sqrt(x.size()) < 1) {
			parentkey.pop_back();
			continue;
		}

		double score = fabs(xTc) + ball.radius*sqrt(x.size());
		if (score >= 1) {
			safeset.push_back(node{parentkey, x, score});
		}

		if (j < d-1 && depth+1 < maxdepth) {
			safeprune(ball, safeset, parentkey, x, depth+1, j+1);
		}

		parentkey.pop_back();
	}
}

void safeprune2(const ball &ball, vector<node> &safeset, int &n_branch_opened){
    double score;
    vector<int> branch_not_screened;
    // Here
    branch_not_screened.reserve(d);
    n_branch_opened = 0;
    for (int i = 0; i < d; i++){
        vector<int> x = Z[i];

		double xTc = 0;
		double p = 0;
		double m = 0;
		for (int k = 0; k < (int)x.size(); k++) {
			double val = ball.center[x[k]];
			xTc += val;
			(val > 0) ? p += val : m += val;
		}

		if (max(p, -m) + ball.radius*sqrt(x.size()) < 1) {
			continue;
		}else{
			n_branch_opened++;
		}
		branch_not_screened.push_back(i);

		score = fabs(xTc) + ball.radius*sqrt(x.size());
		if (score >= 1) {
			safeset.push_back(node{vector<int>{i}, x, score});
		}

        for (int b = 0; b < branch_not_screened.size(); b++){
            int j = branch_not_screened[b];
            if (i==j) continue;
            vector<int> xx; //and
            set_intersection(x.begin(), x.end(), Z[j].begin(), Z[j].end(), back_inserter(xx));
            double xxTc = 0;
            for (int k = 0; k < xx.size(); k++){
                xxTc += ball.center[xx[k]];
            }
            score = fabs(xxTc) + ball.radius*sqrt(xx.size());
            if (score >= 1) {
				safeset.push_back(node{vector<int>{i, j}, xx, score});
            }
        }
    }
    cout << "# branch opened: " << n_branch_opened << endl;
}

// void read(int argc, char **argv) {
// 	int i;
// 	T = 100;
// 	maxdepth = 3;
// 	f = 100;
// 	b = 0;
// 	char filename[1024];
// 	for (i = 1; i < argc; i++) {
// 		if (argv[i][0] != '-') break;
// 		if (++i >= argc)
// 			exit(1);
// 		switch (argv[i-1][1]) {
// 		case 'T':
// 			T = atoi(argv[i]);
// 			break;
// 		case 'D':
// 			maxdepth = atoi(argv[i]);
// 			break;
// 		case 'F':
// 			f = atoi(argv[i]);
// 			break;
// 		case 'B':
// 			b = atoi(argv[i]);
// 			break;
// 		default:
// 			cout << "unknown option" << endl;
// 			exit(1);
// 			break;
// 		}
// 	}

// 	if (i >= argc) exit(1);
// 	strcpy(filename, argv[i]);
// 	read_data(filename);
// }

void read(int argc, char **argv) {
    int i;
    nlambda = 100;
    lambdaMinRatio = 0.01;
    maxSelectedFeatures = 100;
    useBias = 1;
    F = 100;
    eps = 1e-8;
    char filename[1024];
    for (i = 1; i < argc; i++) {
        if (argv[i][0] != '-') break;
        if (++i >= argc) exit(1);
        if (strcmp(argv[i-1], "-nlambda") == 0){
            nlambda = atoi(argv[i]);
        } else if (strcmp(argv[i-1], "-lambdaMinRatio") == 0){
            lambdaMinRatio = atof(argv[i]);
        } else if (strcmp(argv[i-1], "-maxSelectedFeatures") == 0){
             maxSelectedFeatures = atoi(argv[i]);
        }else if (strcmp(argv[i-1], "-useBias") == 0){
            useBias = atoi(argv[i]);
        } else if (strcmp(argv[i-1], "-F") == 0){
            F = atoi(argv[i]);
        } else if (strcmp(argv[i-1], "-eps") == 0){
            eps = atof(argv[i]);
        } else if (strcmp(argv[i-1], "-pathResults") == 0){
        	pathResults = argv[i];
        } else {
            cout << "unknown option" << endl;
            exit(1);
            break;
        }
    }
    if (i >= argc) exit(1);
    strcpy(filename, argv[i]);
    read_data(filename);
}

char* readline(FILE *input) {
	int max_line_len = 1024;
	char* line = (char *)malloc(max_line_len*sizeof(char));

	int len;
	if (fgets(line,max_line_len,input) == NULL) return NULL;

	while (strrchr(line, '\n') == NULL) {
		max_line_len *= 2;
		line = (char *) realloc(line,max_line_len);
		len = (int) strlen(line);
		if (fgets(line+len, max_line_len-len, input) == NULL) break;
	}
	return line;
}

void read_data(char *filename) {
	FILE *fp = fopen(filename, "r");
	if (fp == NULL) {
		fprintf(stderr, "Cannot open input file\n");
		exit(1);
	}

	vector< vector<int> > rowZ;

	n = 0;
	d = 0;
	double Ysum = 0;
	char* line = NULL;
	while ((line = readline(fp)) != NULL) {
		rowZ.push_back(vector<int>());

		char *y = strtok(line, " \t\n");
		Y.push_back(atof(y));
		Ysum += Y[n];

		while(1) {
			char *idx = strtok(NULL, ":");
			char *val = strtok(NULL, " \t");
			if (val == NULL) break;

			int j = atoi(idx);
			rowZ[n].push_back(j);
			if (j > d) d = j;
		}
		n++;
	}
	Ysum /= n;
	for (int i = 0; i < n; i++) {
		Y[i] -= Ysum;
	}

	Z.reserve(d);
	for (int j = 0; j < d; j++) {
		Z.push_back(vector<int>());
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < (int)rowZ[i].size(); j++) {
			int idx = rowZ[i][j]-1;
			Z[idx].push_back(i);
		}
	}
}
