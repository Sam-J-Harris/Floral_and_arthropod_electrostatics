function [EC, maxC, minC, meanC] = FC_err(Vex,Vap,etoly)
% = finds the error between the exact and AAA-LS solutions
EC = min(reftn1(Vex,Vap), reftn2(Vex,Vap,etoly)); % takes minimum of below error functions. EC is a matrix of the error at each point on the meshgrid
maxC = max(EC,[],'all'); minC = min(EC,[],'all'); meanC = mean(EC(~isnan(EC)),'all'); % max, min and mean of errors across the meshgrid
end

%% Appendix: Error functions
% Error function 1: Mean Relative Error
function EC = reftn1(Vex, Vap)
EC = abs(Vex-Vap)./abs(mean(Vex)); % relerr = abs(exact soln - approx soln)/abs(mean(exact soln across the meshgrid))
end

% Error function 2: Relative or Absolute Error
function EC = reftn2(Vex, Vap, etoly)
EC = abs(Vex-Vap)./abs(Vex); % relerr = abs(exact - approx)/abs(exact). Difficulty if either exact=0 or approx = 0.

for k=1:size(EC,1)
    if EC(k)>etoly % if the relerr is above the tolerance (if approx=0 then relerr = 1 > etoly)
        EC(k) = abs(Vex(k)-Vap(k)); % abserr = abs(exact - approx)
    end
end
end
