# Script that trains a classifier to predict the condition on previously unseen runs.
# Trains with different types of classifier, different ROIs, and different left-out runs.
# Appends the job ids to a file jobs.txt for convenience / keeping track of stuff
#

methods=( "cvglmnet" "patternnet" )
masknames=( "hippocampus" "ofc" "striatum" "vmpfc" "bg" "pallidum" "visual" "motor" "sensory" )
runs_to_leaveout=(1 2 3 4 5 6 7 8 9)
z_scores=( "z-none" "z-run" "z-run-voxel" )
predict_what="condition"
trials="[1:24]"

for method in "${methods[@]}"
do
    echo Method $method

    for z_score in "${z_scores[@]}"
    do
        echo Z-score $z_score

        for maskname in "${masknames[@]}"
        do
            mask="masks/$maskname.nii"
            echo Mask $maskname, file = $mask

            for test_run in "${runs_to_leaveout[@]}" 
            do
                training_runs="[1:$((test_run - 1)) $((test_run + 1)):9]"
                echo Leave out run $test_run, train with $training_runs

                outfileprefix="classify_${method}_${z_score}_${maskname}_leaveout_run_${test_run}"
                echo File prefix = $outfileprefix

                # send the job to NCF
                #
                matlab_fn_call="classify('$method', '$mask', '$predict_what', $training_runs, $trials, $test_run, $trials, '$z_score')"
                echo $matlab_fn_call
                sbatch_output=`sbatch -p ncf --mem 25000 -t 20-18:20 -o ${outfileprefix}_%j.out --mail-type=END --wrap="matlab -nodisplay -nosplash -nojvm -r $\"$matlab_fn_call;   exit\""`
                # for local testing
                #sbatch_output=`echo Submitted batch job 88725417`
                echo $sbatch_output

                # Append job id to jobs.txt
                #
                sbatch_output_split=($sbatch_output)
                job_id=${sbatch_output_split[3]}
                echo classify_leaveout_run.sh: $matlab_fn_call -- ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
            done
        done
    done
done

