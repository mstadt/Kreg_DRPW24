% Insulin function
function C_insulin = get_Cinsulin(t_insulin)
      % C_insulin unites are nanomole/L
    % C_insulin units are nanomole/L
    if (t_insulin <= 0)
        %    disp('between meal')
        C_insulin = 22.6/1000;
    elseif (t_insulin > 0) && (t_insulin<(1.5*60))
        %disp('start of meal')
        C_insulin = ((325 - 22.6)/(1.5*60)*(t_insulin) + 22.6)/1000;
    elseif (t_insulin >= (1.5*60)) && t_insulin<(6*60)
        %disp('end of meal')
        C_insulin = ((22.6-325)/((6-1.5)*60)*(t_insulin-6*60) + 22.6)/1000;
    elseif (t_insulin>=(6*60))
        C_insulin = 22.6/1000;
    else
        disp('something weird happened')
        disp(t_insulin)
    end
end % get_C_insulin

