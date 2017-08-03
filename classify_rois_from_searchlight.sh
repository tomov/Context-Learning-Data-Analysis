# convenience script to run the classifier in parallel with different options
#

mkdir output

echo ---------------- >> jobs.txt
echo --- Here we go classify_rois_from_searchlight.sh >> jobs.txt
echo ---------------- >> jobs.txt

methods=( "cvglmnet" )
#masks/glm0_searchlight_sphere_t=10.225_extent=27_roi=Occipital_Inf_R_peak=[36_-84_-16].nii
#masks/glm0_searchlight_sphere_t=9.098_extent=27_roi=Temporal_Inf_L_peak=[-42_-62_-10].nii
#masks/glm0_searchlight_sphere_t=8.656_extent=27_roi=Calcarine_R_peak=[16_-92_2].nii
#masks/glm0_searchlight_sphere_t=8.419_extent=27_roi=Angular_R_peak=[28_-64_46].nii
#masks/glm0_searchlight_sphere_t=6.968_extent=27_roi=Insula_R_peak=[44_22_-4].nii
#masks/glm0_searchlight_sphere_t=6.474_extent=27_roi=Frontal_Mid_2_R_peak=[36_32_22].nii
#masks/glm0_searchlight_sphere_t=5.542_extent=20_roi=Frontal_Inf_Oper_R_peak=[56_18_38].nii
#masks/glm0_searchlight_sphere_t=5.417_extent=27_roi=Frontal_Mid_2_R_peak=[34_4_54].nii
#masks/glm0_searchlight_sphere_t=6.153_extent=27_roi=Precentral_L_peak=[-48_-2_38].nii
#masks/glm0_searchlight_sphere_t=5.144_extent=27_roi=Frontal_Inf_Tri_L_peak=[-42_30_18].nii
#masks/glm0_searchlight_sphere_t=4.666_extent=27_roi=Frontal_Inf_Oper_L_peak=[-50_18_34].nii
masknames=( "glm0_searchlight_sphere_t=8.419_extent=27_roi=Angular_R_peak=[28_-64_46]"  )
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

            outfileprefix="output/classify_${method}_${z_score}_${maskname}_leaveout_run_5"
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

