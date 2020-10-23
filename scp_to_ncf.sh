#rsync -avh --filter='- /*/*/' --include="*/" --include='*.m'  --exclude='*' --rsh="/usr/local/bin/sshpass -p Mainatati13_!_ ssh -o StrictHostKeyChecking=no -l mtomov13" . mtomov13@ncflogin.rc.fas.harvard.edu:/ncf/gershman/Lab/scripts/matlab/ConLearn-ISL/ 
rsync -avh --filter='- /*/*/' --include="*/"  --include='*.sh' --include='*.m'  --exclude='*' --rsh="/usr/local/bin/sshpass -p Mainatati13_!_ ssh -o StrictHostKeyChecking=no -l mtomov13" . mtomov13@ncflogin.rc.fas.harvard.edu:/ncf/gershman/Lab/scripts/matlab/ConLearn-ISL/ 

#scp *.m mtomov13@ncflogin.rc.fas.harvard.edu:/ncf/gershman/Lab/scripts/matlab/ConLearn/
#scp rsa_sanity.m mtomov13@ncflogin.rc.fas.harvard.edu:/ncf/gershman/Lab/scripts/matlab/ConLearn/
#scp mtomov13@ncflogin.rc.fas.harvard.edu:/ncf/gershman/Lab/ConLearn/rsaOutput/rsa1/searchlight*.mat ../rsaOutput/rsa1/
#scp mtomov13@ncflogin.rc.fas.harvard.edu:/ncf/gershman/Lab/scripts/matlab/ConLearn/rdms/searchlight*.mat rdms/M1M2M1_4mm_ordered/
#scp mtomov13@ncflogin.rc.fas.harvard.edu:/ncf/gershman/Lab/scripts/matlab/ConLearn/slurm_matlab.sh .
#scp mtomov13@ncflogin.rc.fas.harvard.edu:/ncf/gershman/Lab/scripts/matlab/ConLearn/rdms_searchlight.sh .
#scp mtomov13@ncflogin.rc.fas.harvard.edu:/ncf/gershman/Lab/ConLearn/glmOutput/mean.nii ../glmOutput/
#scp results/fit_params_results_* mtomov13@ncflogin.rc.fas.harvard.edu:/net/rcss11/srv/export/ncf_gershman/share_root/Lab/scripts/matlab/ConLearn/results/
#scp atlases/* mtomov13@ncflogin.rc.fas.harvard.edu:/net/rcss11/srv/export/ncf_gershman/share_root/Lab/scripts/matlab/ConLearn/atlases/
#scp results/fit_params_results_*collins* mtomov13@ncflogin.rc.fas.harvard.edu:/net/rcss11/srv/export/ncf_gershman/share_root/Lab/scripts/matlab/ConLearn/results/
#scp results/fit_params_results_rev*.mat mtomov13@ncflogin.rc.fas.harvard.edu:/net/rcss11/srv/export/ncf_gershman/share_root/Lab/scripts/matlab/ConLearn/results/
#scp masks/*.nii mtomov13@ncflogin.rc.fas.harvard.edu:/net/rcss11/srv/export/ncf_gershman/share_root/Lab/scripts/matlab/ConLearn/masks/
