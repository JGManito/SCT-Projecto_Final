function [satellite_positions_minimized] = minimize_PDOP(satellite_positions,position_estimate,n_input)

NSats = size(satellite_positions,1);



list_sats = 1:NSats;

n=n_input;%Number of satellites to use as a sub-constellation
    %Clear variables
    
    H=[];
    combs_DOP_SVN=[];
    PDOP=[];
    HDOP=[];
    
    %List all possible n satellite combinations
    combs_DOP=nchoosek(list_sats,n);
    
    %Iterate all the permutations and get the H matrix for each one of them
    for j=1:size(combs_DOP,1)
        positions_aux=[];
        %Get the positions of the satellite subset
        positions_aux=[positions_aux;satellite_positions(combs_DOP(j,:),:)];
        
        %Get H matrix for the permutation
        for i=1:n
            
            D_estimate(i) = norm(positions_aux(i,:)-position_estimate(1:3));
            
            G(i,1) = -(positions_aux(i,1)-position_estimate(1))/D_estimate(i);
            G(i,2) = -(positions_aux(i,2)-position_estimate(2))/D_estimate(i);
            G(i,3) = -(positions_aux(i,3)-position_estimate(3))/D_estimate(i);
            G(i,4) = 1;
            
        end

        H = inv(G'*G);
        
        %Get DOP matrix
        DOP_matrix=inv(transpose(H)*H);
        
        PDOP_current = sqrt(trace(DOP_matrix(1:3,1:3)));
        
        %Compute PDOP from DOP matrix
        if PDOP_current == 0
            PDOP_current = 9999;
        elseif isreal(PDOP_current)
            PDOP(j)=PDOP_current;
        else
            PDOP(j) = 9999;
        end
    end
    
    %Get minimum PDOP and HDOP
    [min_PDOP,min_PDOP_index]=min(PDOP);

combs_DOP=nchoosek(list_sats,n);

%Get the pseudoranges and positions of the satellite subset
satellite_positions_minimized=satellite_positions(combs_DOP(min_PDOP_index,:),:);


end
