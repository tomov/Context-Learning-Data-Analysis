# convenience script to run the classifier in parallel with different options
#

mkdir output

echo ---------------- >> jobs.txt
echo --- Here we go classify.sh >> jobs.txt
echo ---------------- >> jobs.txt

methods=( "cvglmnet" )
masknames=( "glm0 searchlight sphere t=10.225 extent=27 roi=Occipital_Inf_R peak=[36 -84 -16]"  "glm0 searchlight sphere t=9.098 extent=27 roi=Temporal_Inf_L peak=[-42 -62 -10]"  "glm0 searchlight sphere t=8.656 extent=27 roi=Calcarine_R peak=[16 -92 2]"  "glm0 searchlight sphere t=8.419 extent=27 roi=Angular_R peak=[28 -64 46]"  "glm0 searchlight sphere t=6.968 extent=27 roi=Insula_R peak=[44 22 -4]"  "glm0 searchlight sphere t=6.474 extent=27 roi=Frontal_Mid_2_R peak=[36 32 22]"  "glm0 searchlight sphere t=5.542 extent=20 roi=Frontal_Inf_Oper_R peak=[56 18 38]"  "glm0 searchlight sphere t=5.417 extent=27 roi=Frontal_Mid_2_R peak=[34 4 54]"  "glm0 searchlight sphere t=6.153 extent=27 roi=Precentral_L peak=[-48 -2 38]"  "glm0 searchlight sphere t=5.144 extent=27 roi=Frontal_Inf_Tri_L peak=[-42 30 18]"  "glm0 searchlight sphere t=4.666 extent=27 roi=Frontal_Inf_Oper_L peak=[-50 18 34]" )
#masknames=( "hippocampus" "ofc" "striatum" "vmpfc" "bg" "pallidum" "visual" "v1" "motor" "sensory" )
z_scores=( "z-none" )
predict_what="condition"
training_runs="[1:4 6:9]"
training_trials="[1:20]"
test_runs="[5]"
test_trials="[1:20]"

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
            #sbatch_output=`sbatch -p ncf --mem 25000 -t 20-18:20 -o ${outfileprefix}_%j.out --mail-type=END --wrap="matlab -nodisplay -nosplash -nojvm -r $\"$matlab_fn_call;   exit\""`
            # for local testing
            sbatch_output=`echo Submitted batch job 88725417`
            echo $sbatch_output

            # Append job id to jobs.txt
            #
            sbatch_output_split=($sbatch_output)
            job_id=${sbatch_output_split[3]}
            echo classify.sh: $matlab_fn_call -- ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
        done
    done
done

