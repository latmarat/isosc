function [c,conv] = bisect(a,b,s1,s2,f,m,itermax,tol,xFlag)
    % find solution for isoSC using bisection method 
    
    ii = 0;
    conv = 0;
    
    % start the loop
    while ii <= itermax
        c = (a+b)/2;
        
        % solve i.t.o x1 or x2
        if xFlag == 1
            err = func2solve1(s1,s2,f,m,c);
        elseif xFlag == 2
            err = func2solve2(s1,s2,f,m,c);
        end
        
        % convergence
        if abs(err) <= tol
            conv = 1;
            break
        end
        
        % otherwise find the sign
        if xFlag == 1
            errSign = sign(func2solve1(s1,s2,f,m,a));
        elseif xFlag == 2
            errSign = sign(func2solve2(s1,s2,f,m,a));
        end
  
        % set new interval bounds
        if sign(err) == errSign
            a = c;
        else
            b = c;
        end
        ii = ii + 1;
    end
    fprintf('Bisect has done %d iterations\n',ii)
    fprintf('Final error: %e\n',err)
end

function err = func2solve1(s1,s2,f1,m,x1)
    ac = 5/3;
    x2 = (1.0-f1*x1)/(1.0-f1);
    err = s1*x1^m/(ac-x1) - s2*x2^m/(ac-x2);
end

function err = func2solve2(s1,s2,f2,m,x2)
    ac = 5/3;
    x1 = (1.0-f2*x2)/(1.0-f2);
    err = s1*x1^m/(ac-x1) - s2*x2^m/(ac-x2);
end
