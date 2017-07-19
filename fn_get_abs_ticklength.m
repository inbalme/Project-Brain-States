% setting TickLength in absolute value
function ticklen=fn_get_abs_ticklength(hax, abslen)
%hax is the axis handle
%abslen is the desired target length (in figure units)

rect = get(hax,'Position');     % Axis position (relative to figure)
hfig = get(hax,'Parent');       % Handle to parent figure
rectfig = get(hfig, 'Position'); % Figure position (in pixels or other figure units)
width = rect(3) * rectfig(3);    % Absolute width of axis (in pixels or other figure units)
height = rect(4) * rectfig(4);   % Absolute height of axis (in pixels or other figure units)
    axislen = max([height,width]);    % Get longest axis
    ticklen = abslen/axislen;         % Fix length
    set(hax,'TickDir','out','TickLength',ticklen*[1 1]);
end