function rdms_searchlight_batch(start_idx, end_idx, r)

batch_size = 1000;
N = ceil((end_idx - start_idx + 1) / batch_size);

for i = start_idx:batch_size:end_idx
    s = i;
    e = i + batch_size - 1;
    fprintf('RANGE %d - %d\n', s, e);
    rdms_searchlight(s, e, r);
end
