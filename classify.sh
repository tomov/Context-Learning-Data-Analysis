# convenience script to run the classifier in parallel with different options
#

mkdir output

echo ---------------- >> jobs.txt
echo --- Here we go classify.sh >> jobs.txt
echo ---------------- >> jobs.txt

methods=( "cvglmnet" )
masknames=( "hippocampus" )
#masknames=( "hippocampus" "ofc" "striatum" "vmpfc" "bg" "pallidum" "visual" "v1" "motor" "sensory" )
z_scores=( "z-none" )
predict_what="condition"
training_runs="[1:9]"
training_trials="[1:16]"
test_runs="[1:9]"
test_trials="[17:20]"

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

            outfileprefix="output/classify_${method}_${z_score}_${maskname}_leaveout_trials_17-20"
            echo File prefix = $outfileprefix

            # send the job to NCF
            #
            matlab_fn_call="classify('$method', '$mask', '$predict_what', $training_runs, $training_trials, $test_runs, $test_trials, '$z_score')"
            echo $matlab_fn_call
            sbatch_output=`sbatch -p ncf --mem 25000 -t 20-18:20 -o ${outfileprefix}_%j.out --mail-type=END --wrap="matlab -nodisplay -nosplash -nojvm -r $\"$matlab_fn_call;   exit\""`
            # for local testing
            #sbatch_output=`echo Submitted batch job 88725417`
            echo $sbatch_output

            # Append job id to jobs.txt
            #
            sbatch_output_split=($sbatch_output)
            job_id=${sbatch_output_split[3]}
            echo classify.sh: $matlab_fn_call -- ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
        done
    done
done

