function [] = forces_drawer(node, force, magn, colour)
    %%%
    %%% Plot force (attenzione alla forma dei vettori di input!)
    %%%
    
    %% Horizontal forces
    i = 2;
    if force(i) ~= 0
        p1 = [(node(force(1),2)-(force(i)/(2000000*magn))) node(force(1),3)];
        p2 = [node(force(1),2) node(force(1),3)];
        dp = p2-p1;
        quiver(p1(1),p1(2),dp(1),dp(2),0, 'Color', colour)
        text((p1(1)+p2(1)-magn)/2,(p1(2)+p2(2)+magn)/2, sprintf('%.0f kN',force(i)/1000)) 
    end
    clear i    
    
    %% Vertical forces
    i = 3;
    if force(i) ~= 0
        p1 = [node(force(1),2) (node(force(1),3)-(force(i)/(2000000*magn)))];
        p2 = [node(force(1),2) node(force(1),3)];
        dp = p2-p1;
        quiver(p1(1),p1(2),dp(1),dp(2),0, 'Color', colour)
        text((p1(1)+p2(1))/2,(p1(2)+p2(2)+magn)/2, sprintf('%.0f kN',abs(force(i))/1000)) 
    end 
    clear i
  %https://it.mathworks.com/matlabcentral/answers/160487-how-can-i-draw-a-line-with-arrow-head-between-2-data-points-in-a-plot
end