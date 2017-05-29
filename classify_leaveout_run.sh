# Script that trains a classifier to predict the condition on previously unseen runs.
# Trains with different types of classifier, different ROIs, and different left-out runs.
#

methods=( "cvglmnet" "patternnet" )
masknames=( "hippocampus" "ofc" "striatum" "vmpfc" "bg" "pallidum" "visual" "motor" "sensory" )
runs_to_leaveout=(1 2 3 4 5 6 7 8 9)
predict_what="condition"
trials="[1:24]"

for method in "${methods[@]}"
do
    echo Method $method
    for maskname in "${masknames[@]}"
    do
        mask="masks/$maskname.nii"
        echo Mask $maskname, file = $mask

        for test_run in "${runs_to_leaveout[@]}" 
        do
            training_runs="[1:$((run - 1)) $((run + 1)):9]"
            echo Leave out run $test_run, train with $training_runs

            outfileprefix="classify_${method}_${maskname}_leaveout_run_${test_run}"
            echo File prefix = $outfileprefix

            sbatch -p ncf --mem 25000 -t 20-18:20 -o ${outfileprefix}_%j.out --mail-type=END --wrap="matlab -nodisplay -nosplash -nojvm -r $\"classify('$mask', '$predict_what', $trials, $training_runs, $trials, $test_run);   exit\""
        done
    done
done

