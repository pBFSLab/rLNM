function cfg = specificity_config()
% Shared settings for the specificity pipeline.
%
% SIDED is used by S5_permutation_fdr and specificity_rethreshold, so both
% always report the same side. Edit it here only.

cfg.fdr_level = 0.05;
cfg.sided = 'one';            % 'one' or 'two'
cfg.n_permutations = 10000;
cfg.n_workers = 8;
end
