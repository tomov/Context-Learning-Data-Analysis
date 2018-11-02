rsa_behav_sanity;

rsa_neural_sanity;


[Rho, H, T, P, all_subject_rhos] = ccnl_match_rdms(N, B, c1);

[table_Rho, table_H, table_T, table_P, all_subject_rhos] = rdms_second_order(metadata, S, M, c2, false, [], []);


assert(isequal(Rho(:,1), table_Rho(:,1)));
