% This function generates a character string using the current date in the
% following format: yymmdd

function t = timestamp(varargin)
    c = clock; % current date and time array
    year = num2str(c(1));
    month = num2str(c(2));
    if length(num2str(c(2))) == 1
        month = strcat('0',month);
    end
    day = num2str(c(3));
    if length(day) == 1
        day = strcat('0',day);
    end
    t = strcat(year(3:4),month,day);
end