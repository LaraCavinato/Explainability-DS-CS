function [ci] = concordance_index(targets, risks, censored, ties)
% targets: the respone, the 'time to event'
% risks: the risks that my risk model have assigned to that patient. 
%        A patient with higher risk score should have a shorter time to 
%        event censored is a boolean that says if the relative time to 
%        event is observed or not. 1 for data that are observed and 0 for 
%        data censored.
% ties: flag used to choose how to handle the ties case, in particular
%       ties = 1 means that you consider the ties case and if 0 you ignore them. 
  [targets,i] = sort(targets,'descend');
  risks = risks(i);
  censored = censored(i);
  n = length(targets);
  total  = 0;
  norm_z = 0;
  for j=1:n
      for k=(j+1):n
          if targets(j) ~= targets(k) 
              if(censored(k)==1)
                h = step_function(risks(j) - risks(k),ties);
                total = total + h;
                norm_z = norm_z + 1;
              end
          end
      end
  end
  ci = total / norm_z;
end
 function h = step_function(diff, ties)
    if diff < 0
        h = 1;
    elseif ((diff == 0) & (ties==1))
        h = 0.5;
    else
        h = 0;
    end
 end