function [] = restraints_drawer(x,y,fact_x,fact_y,type,colour)
    %%%
    %%% 1 for rollers
    %%% 2 for hinges
    %%% joints not jet implemented
    %%%

    % Points coordinates
    if type == 1    % roller case
        a = [x-fact_x y-(fact_y*1.5)];
        b = [x-(fact_x/2) y-fact_y];
        c = [x y];
        d = [x+(fact_x/2) y-fact_y];
        e = [x+fact_x y-(fact_y*1.5)];
    elseif type == 2    % hinge case
        a = [x-fact_x y-fact_y];
        b = [x-(fact_x/2) y-fact_y];
        c = [x y];
        d = [x+(fact_x/2) y-fact_y];
        e = [x+fact_x y-fact_y];
    elseif type == 3 % horizontal joint case
    elseif type == 4 % vertical joint case
    end
    
    
    % Lines coordinates
    xs = [b(1) c(1) d(1) b(1) a(1) e(1)];
    ys = [b(2) c(2) d(2) b(2) a(2) e(2)];
    
    
    % Point plot
    plot(xs(1:4),ys(1:4), 'Color', colour, 'HandleVisibility', 'off')
    plot(xs(5:6),ys(5:6), 'Color', colour, 'HandleVisibility', 'off')
    
    %Field plot
    y_field = [a(2) a(2)-(fact_y/2)];
    field_lenght = (xs(6) - xs(5))/4;
    for i = 0:4
        x_field = [xs(5)+field_lenght*i xs(5)+field_lenght*(i-1)];
        plot(x_field, y_field,'Color', colour, 'HandleVisibility', 'off')
    end 
    clear i
    

    % Roller wheels
    if type == 1
        r = (b(2)-a(2))/2;
        y = b(2) - r;
        x1 =(b(1) + (d(1) - b(1))/4);
        x2 =(d(1) - (d(1) - b(1))/4);
        circle(x1,y,r,colour)
        circle(x2,y,r,colour)
    end

    function circle(x,y,r,colour)
        circ=0:0.01:2*pi; 
        xc=r*cos(circ);
        yc=r*sin(circ);
        plot(x+xc,y+yc, 'Color', colour, 'HandleVisibility','off');
            
    end
end
