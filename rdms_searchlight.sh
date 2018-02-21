# run searchlight for ranges of voxels
#

mkdir output

echo ---------------- >> jobs.txt
echo --- Here we go rdms_searchlight >> jobs.txt
echo ---------------- >> jobs.txt

batch_size=1000;
r=1.814;

for batch in {1..10}
do
	start_idx=$(((batch - 1) * batch_size + 1))
	end_idx=$((batch * batch_size))

	outfileprefix="output/rdms_searchlight_${start_idx}_${end_idx}_${r}"
	echo File prefix = $outfileprefix

	rdms_searchlight_call="rdms_searchlight(${start_idx},${end_idx},${r})"
	echo $rdms_searchlight_call

    sbatch_output=`CMD="$rdms_searchlight_call" sbatch -p ncf --mem 25000 -t 20-18:20 -o ${outfileprefix}_%j.out --mail-type=END slurm_matlab.sh`
	#sbatch_output=`sbatch -p ncf --mem 25000 -t 20-18:20 -o ${outfileprefix}_%j.out --mail-type=END --wrap="matlab -nodisplay -nosplash -nojvm -r $\"$rdms_searchlight_call; exit\""`
        #sbatch_output=`echo Submitted batch job 88725418`
	echo $sbatch_output

	# Append job id to jobs.txt
	#
	sbatch_output_split=($sbatch_output)
	job_id=${sbatch_output_split[3]}
 	echo rdms_searchlight.sh: ${outfileprefix}_${job_id}.out -- $sbatch_output >> jobs.txt
done
