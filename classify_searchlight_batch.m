function classify_searchlight_batch(start_idx, end_idx, r, event)
% Call classify_searchlight in batches so it saves every now and then
% similar to rdms_searchlight_batch
%

batch_size = 10000;
N = ceil((end_idx - start_idx + 1) / batch_size);

for i = start_idx:batch_size:end_idx
    s = i;
    e = i + batch_size - 1;
    fprintf('RANGE %d - %d\n', s, e);
    classify_searchlight(s, e, r, event);
end
