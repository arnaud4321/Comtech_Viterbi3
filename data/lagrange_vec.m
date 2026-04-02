function sig_out = lagrange_vec(sig_in, t_target, order)
    % sig_in   : Vecteur du signal original (1xN)
    % t_target : Instants désirés (ex: 1, 1.66, 2.33...)
    % order    : Ordre du polynôme (3 = quadratique, 4 = cubique)
    
    n_out = length(t_target);
    sig_out = zeros(size(t_target));
    half_order = floor(order / 2);
    
    for i = 1:n_out
        t = t_target(i);
        k = floor(t); % Indice de base
        
        % Extraction du voisinage (support)
        % On s'assure de ne pas déborder du vecteur original
        idx_start = max(1, k - half_order + 1);
        idx_end   = min(length(sig_in), idx_start + order - 1);
        idx_start = max(1, idx_end - order + 1); % Recalage si bord droit
        
        x_data = idx_start:idx_end;
        y_data = sig_in(x_data);
        
        % Calcul de Lagrange pour le point t
        val = 0;
        for j = 1:length(x_data)
            L = 1;
            for m = 1:length(x_data)
                if j ~= m
                    L = L * (t - x_data(m)) / (x_data(j) - x_data(m));
                end
            end
            val = val + y_data(j) * L;
        end
        sig_out(i) = val;
    end
end