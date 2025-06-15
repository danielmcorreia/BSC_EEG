function [pvals_, pvals_adjust_] = get_pvals(chance, real)

    % get pvalues
    pvals=stat_surrogate_pvals(chance, real','one');

    % Adjust for multiple comparisons
    p1=pval_adjust(pvals(:,1), 'fdr'); p2=pval_adjust(pvals(:,2), 'fdr'); p3=pval_adjust(pvals(:,3), 'fdr');
    pvals_adjust_(:,1) = p1; pvals_adjust_(:,2) = p2; pvals_adjust_(:,3) = p3;

    % or not
    p1=pvals(:,1); p2=pvals(:,2); p3=pvals(:,3);
    pvals_(:,1) = p1; pvals_(:,2) = p2; pvals_(:,3) = p3;