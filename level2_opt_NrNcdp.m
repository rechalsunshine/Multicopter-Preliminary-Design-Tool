%%% This function is a wrap of the level2_opt function, it outputs the MOTM

function MTOM=level2_opt_NrNcdp(DL)
global  result t count FUNoutput
err_per=level2_NrNcdp(DL);
% if abs(FUNoutput(19))<0.05
    MTOM=FUNoutput(18);
result(t,:)=FUNoutput;
t=t+1;
err_100=FUNoutput(19)*100;
fprintf('current count: %d, MTOM: %f g, error: %.1f%% \n', count, MTOM, err_100);
% else
%         MTOM=100000;
%         fprintf('current count: %d, invalid results...\n', count);
% end
end