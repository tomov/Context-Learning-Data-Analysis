# Script that trains a classifier to predict the condition from trial onset betas
# on all trials on all runs, for a single subject.
# Train a separate classifier with different types of classifiers, different ROIs, and different subjects.
# Appends the job ids to a file jobs.txt for convenience / keeping track of stuff
#

echo ---------------- >> jobs.txt
echo --- Here we go classify_each_subject >> jobs.txt
echo ---------------- >> jobs.txt

methods=( "cvglmnet" "patternnet" )
masknames=( "hippocampus" "ofc" "striatum" "vmpfc" "bg" "pallidum" "visual" "motor" "sensory" )
z_scores=( "z-none" "z-run" "z-run-voxel" )
goodSubjects=( 1 2 4 6 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 )  # getGoodSubjects()
predict_what="condition"
runs="[1:9]"
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

            for subj in "${goodSubjects[@]}"
            do
                echo Subject $subj

                outfileprefix="classify_train_${method}_${z_score}_${maskname}_subj_${subj}"
                echo File prefix = $outfileprefix

                # send the job to NCF
                #
                matlab_fn_call="classify_train('$method', $runs, $trials, [$subj], '$mask', '$predict_what', '$z_score')"
                echo $matlab_fn_call
                sbatch_output=`sbatch -p ncf --mem 25000 -t 0-18:20 -o ${outfileprefix}_%j.out --mail-type=END --wrap="matlab -nodisplay -nosplash -nojvm -r $\"$matlab_fn_call;   exit\""`
                # for local testing
                #sbatch_output=`echo Submitted batch job 88725418`
                echo $sbatch_output

                # Append job id to jobs.txt
                #
                sbatch_output_split=($sbatch_output)
                job_id=${sbatch_output_split[3]}
                echo classify_each_subject.sh: $matlab_fn_call -- ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
            done
        done
    done
done

