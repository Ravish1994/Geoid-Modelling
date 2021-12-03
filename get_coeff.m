function [C S] = get_coeff(mat,l,m)
for i = 1:length(mat)
    if mat(i,1) == l
        for j = i:i+l
            if mat(j,2) == m
                C = mat(j,3);
                S = mat(j,4);
                break;
            end
        end
        break;
    end
end
end