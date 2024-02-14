function pos = generate_NPs(rSurf, dmin, num_points, pos)
    rho = rand(1,num_points) + rand(1,num_points);
    rho(rho>1)=2-rho(rho>1);
    theta = 2*pi*rand(1,num_points);

    new_pos = [rSurf*rho.*cos(theta); rSurf*rho.*sin(theta)];
    
    new_pos = single(new_pos);
    pos = single(pos);
    
    checkoverlap = true;
    
    while checkoverlap
        dist_between_nps = dist(new_pos', new_pos);
        dist_between_nps = dist_between_nps - diag(diag(dist_between_nps));
        dist_between_nps = dist_between_nps + diag(5000*ones(1,size(new_pos,2)));
        
        overlaps = dist_between_nps < dmin;
        [I, J] = find(overlaps);
        
        if ~isempty(I)
                rho = rand(1,length(I)) + rand(1,length(I));
                rho(rho>1)=2-rho(rho>1);
                theta = 2*pi*rand(1,length(I));
                new_pos(:,I) = [rSurf*rho.*cos(theta); rSurf*rho.*sin(theta)];
        else
            checkoverlap = false;
        end
        
    end
    
    dist_between_nps = dist(new_pos', pos);             % new_pos : 2xN ----pos : 2xM ----- dist : NxM
    dist_between_nps = real(dist_between_nps);

    overlaps = dist_between_nps < dmin;
    [I, J] = find(overlaps);
    new_pos(:,I)=[];
    
    if any(min(dist_between_nps) < dmin)
        disp(find(dist_between_nps(i,:) == min(dist_between_nps(i,:))));
        error('NP overlap')
    end
    
    pos = [pos, new_pos];
    
end