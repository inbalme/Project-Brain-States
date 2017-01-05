%Change capsize of errorbars
%asymmetric=1 - change only size of lower errorbars, asymmetric=2 - change only upper ones,
%asymmetric=0 for both sides
%remove = 1 if want to remove the cap completely.
function fn_errorbar_capsize(errbar_h,capincrease,asymmetric,remove)
% hb = get(errbar_h,'children'); 
% Xdata = get(hb(2),'Xdata');
Xdata = get(errbar_h,'Xdata');
temp = 4:3:length(Xdata);
temp(3:3:end) = [];
if asymmetric==1;
temp(1:2:end)=[];
end
if asymmetric==2;
temp(2:2:end)=[];
end
% xleft and xright contain the indices of the left and right
%  endpoints of the horizontal lines
xleft = temp; xright = temp+1;
% Increase line length by 0.2 units
Xdata(xleft) = Xdata(xleft) - capincrease/2;
Xdata(xright) = Xdata(xright) +capincrease/2;
set(errbar_h,'Xdata',Xdata)
% set(hb(2),'Xdata',Xdata)
end