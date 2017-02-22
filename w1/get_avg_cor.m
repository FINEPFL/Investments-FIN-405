function avg_cor = get_avg_cor(return_table)
    cor_mat = corr(return_table);
    avg_cor = 1/(7 * (7-1)) * (sum(sum(cor_mat)) - sum(diag(cor_mat)));
    
    
