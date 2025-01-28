function [f107,f107p,f107a,Ap] = activity(date,opts)
arguments
    date (1,1) datetime
    opts.range (1,1) int32 {mustBePositive} = 81
end

num_header_lines = 41;
delta_days = round((opts.range - 1) / 2);
id0 = days(date - datetime(1932,1,1)) + num_header_lines - delta_days;

url = 'https://www-app3.gfz-potsdam.de/kp_index/Kp_ap_Ap_SN_F107_since_1932.txt';
lines = strsplit(webread(url),'\n');

vals = nan(2,opts.range);
for id = id0:id0 + opts.range - 1
    data = strsplit(lines{id});
    vals(1,id-id0+1) = str2double(data{24});
    vals(2,id-id0+1) = str2double(data{26});
end

f107 = vals(2, delta_days + 1);
f107p = vals(2, delta_days);
f107a = mean(vals(2, :));
Ap = vals(1, delta_days + 1);
