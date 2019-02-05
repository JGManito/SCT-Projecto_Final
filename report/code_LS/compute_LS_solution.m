function [position_estimate] = compute_LS_solution(satellite_positions,initial_estimate,pseudoranges)

NSat = size(satellite_positions,1);

position_estimate = [initial_estimate(1:3) 0];
previous_estimate = initial_estimate(1:3);

stop = 1;

while stop
    
    for i=1:NSat
        
        D_estimate(i) = norm(satellite_positions(i,:)-position_estimate(1:3));
        
        dPseudorange(i) = pseudoranges(i) - (D_estimate(i)+position_estimate(4));
        
        G(i,1) = -(satellite_positions(i,1)-position_estimate(1))/D_estimate(i);
        G(i,2) = -(satellite_positions(i,2)-position_estimate(2))/D_estimate(i);
        G(i,3) = -(satellite_positions(i,3)-position_estimate(3))/D_estimate(i);
        G(i,4) = 1;
        
    end
    
    dX = inv(G'*G)*G'*dPseudorange';
    
    position_estimate = position_estimate + dX';
        

    if norm(abs(position_estimate(1:3) - previous_estimate)) < 0.001
        stop = 0;
    end
    
    previous_estimate = position_estimate(1:3);



end

