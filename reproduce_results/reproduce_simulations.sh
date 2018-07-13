# Simulations
# -----------

# Generate the simulated data
python ./get_data_preprocessed/simulate_lasso.py

# Create the folders to hold the results
mkdir -p ./results/WHInter/
mkdir -p ./results/zetaIL/
mkdir -p ./results/SPP/
mkdir -p ./results/Blitz/

# Run
for i in 1000,1000 1000,3000 1000,10000 300,1000 10000,1000
do
    IFS=","
    set -- $i
    dat=./data_preprocessed/Bernoulli_n${1}_p${2}_qunif_coefnormal_rs0_nnzd.tsv
    ./../src/train_WHInter -nlambda 100 -lambdaMinRatio 0.01 -maxSelectedFeatures 150 -useBias 1 -useMyMips 0 -typeBound 2 -F 50 -pathResults ./results/WHInter/ ${dat}
    ./../src/train_WHInter -nlambda 100 -lambdaMinRatio 0.01 -maxSelectedFeatures 150 -useBias 1 -useMyMips 1 -typeBound 2 -F 50 -pathResults ./results/WHInter/ ${dat}
    ./../src/train_WHInter -nlambda 100 -lambdaMinRatio 0.01 -maxSelectedFeatures 150 -useBias 1 -useMyMips 2 -typeBound 2 -F 50 -pathResults ./results/WHInter/ ${dat}
    ./../src/train_WHInter -nlambda 100 -lambdaMinRatio 0.01 -maxSelectedFeatures 150 -useBias 1 -useMyMips 2 -typeBound 1 -F 50 -pathResults ./results/WHInter/ ${dat}
    ./../src/train_WHInter -nlambda 100 -lambdaMinRatio 0.01 -maxSelectedFeatures 150 -useBias 1 -useMyMips 2 -typeBound 0 -F 50 -pathResults ./results/WHInter/ ${dat}
    ./src/zetaIL/train_zetaIL -nlambda 100 -lambdaMinRatio 0.01 -maxSelectedFeatures 150 -useBias 1 -useMyMips 0 -F 50 -pathResults ./results/zetaIL/ ${dat}
    ./src/zetaIL/train_zetaIL -nlambda 100 -lambdaMinRatio 0.01 -maxSelectedFeatures 150 -useBias 1 -useMyMips 1 -F 50 -pathResults ./results/zetaIL/ ${dat}
    ./src/zetaIL/train_zetaIL -nlambda 100 -lambdaMinRatio 0.01 -maxSelectedFeatures 150 -useBias 1 -useMyMips 2 -F 50 -pathResults ./results/zetaIL/ ${dat}
    ./src/SPP/train_SPPbreadthFirst -nlambda 100 -lambdaMinRatio 0.01 -maxSelectedFeatures 150 -useBias 1 -F 50 -pathResults ./results/SPP/ ${dat}
    if [ ${2} -ne 10000 ]
    then
        python ./src/Blitz/runBLITZ.py --results_dir ./results/Blitz --data_file ${dat} --nlambda 100 --lambdaMinRatio 0.01 --maxSelectedFeatures 150 --useBias 1 --tol 1e-8
    fi
done

# Plot
Rscript ./analyze_results/analyze_sim.R
Rscript ./analyze_results/analyze_support_recovery.R