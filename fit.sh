outfileprefix="output/fit_params_"
echo File prefix = $outfileprefix

#command="fit_params(0, 1, {'simple_collins'}, 25, 'results/fit_params_results_simple_collins_25nstarts_0-10alpha.mat')"
#command="fit_params(0, 1, {[1 1 1 0 0]}, 25, 'results/fit_params_results_M1M2M3_25nstarts.mat')"
#command="fit_params(0, 1, {[1 1 0 1 0]}, 25, 'results/fit_params_results_M1M2M1_25nstarts.mat')"
#command="fit_params(0, 1, {[1 1 1 0 0]}, 25, 'results/fit_params_results_M1M2M3_25nstarts_tau.mat', 3)"
#command="fit_params(0, 1, {[1 1 0 1 0]}, 25, 'results/fit_params_results_M1M2M1_25nstarts_tau.mat', 3)"
#command="fit_params(0, 1, {[1 1 1 0 0]}, 25, 'results/fit_params_results_M1M2M3_25nstarts_tau_w0.mat', 4)"
#command="fit_params(0, 1, {[1 1 0 1 0]}, 25, 'results/fit_params_results_M1M2M1_25nstarts_tau_w0.mat', 4)"
command="fit_params(0, 1, {'simple_collins'}, 25, 'results/fit_params_results_simple_collins_25nstarts_0-10alpha_Q0.mat', 4)"
echo $command

sbatch_output=`CMD="$command" sbatch -p ncf --mem 50000 -t 0-3:20 -o ${outfileprefix}_%j.out --mail-type=END slurm_matlab.sh`
    #sbatch_output=`echo Submitted batch job 88725418`
echo $sbatch_output

# Append job id to jobs.txt
#
sbatch_output_split=($sbatch_output)
job_id=${sbatch_output_split[3]}
echo fit.sh: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt

