function mi = mutual_info(X,Y,norm)

    H_X = shannon_entropy(X);
    H_Y = shannon_entropy(Y);
    H_XY = joint_entropy(X,Y);

    mi = H_X + H_Y - H_XY;

    if norm
        mi = mi ./ (H_X + H_Y);
    end
     
    mi = H_X;
end

