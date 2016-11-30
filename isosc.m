function [ convFlag, x1, x2, sbar2, check ] = isosc( f1, eta, m  )
    % Implementation of Stringfellow-Parks model
    % constants
    a = 5/3;
    b = 2/3;
    
    % set upper limit for x1 to avoid singularity
    x1lim = fix(100*a)/100.0;

    % tolerance
    tiny = 1.0e-12; % realmin('double');
    
    % hardness and volume fraction
    s2 = 1.0;
    s1 = s2/eta;
    f2 = 1.0 - f1;

    flag2 = 0;
    
    % solve for x1
    [x1,flag1] = bisect(1.0,x1lim,s1,s2,f1,m,1000,tiny,1);
    x2 = (1-f1*x1)/f2;
    
    % if not converget, solve for x2
    if flag1 == 0
        [x2,flag2] = bisect(0.0,1.0,s1,s2,f2,m,1000,tiny,2);
        x1 = (1-f2*x2)/f1;
    end
    
    % find composite hardness ito of x1 and x2
    sbar1 = b*s1*x1^m/(a-x1);
    sbar2 = b*s2*x2^m/(a-x2);
    
    % check if they are the same
    check = sbar2 - sbar1;
%     fprintf('x1: %0.4f,\tx2:%0.4f\n\n',x1,x2)
%     fprintf('sbar1: %0.2f,\tsbar2: %0.2f\n\n',sbar1,sbar2)
%     fprintf('f1x1 + f2x2 = %0.3f\n',f1*x1+f2*x2)

    % convergence flag
    convFlag = sum([flag1,flag2]);
%     if sum([flag1,flag2]) == 0
%         fprintf('Warning! Not converged\n')
%     end

end