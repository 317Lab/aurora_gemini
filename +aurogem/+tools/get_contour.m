function [cx,cy] = get_contour(M,opts)
arguments
    M (2,:) double {mustBeNonempty}
    opts.rank (1,1) int8 {mustBePositive} = 1
end
%#ok<*AGROW>
rank = opts.rank;

lt = size(M,2);
i = 1;
ci_prev = -inf*ones(2,1);
ranks = [];
while true
    li = M(2,i);
    ci = M(:,i+1:i+li);
    yi = mean(ci(2,:));
    yi_prev = mean(ci_prev(2,:));
    if all(ci(:,1)==ci(:,li)) % closed contour
        break
    elseif yi > yi_prev
        ranks = [i,ranks];
    else
        ranks = [ranks,i];
    end
    ci_prev = ci;
    i = i + li + 1;
    if i > lt
        break
    else
    end
end
if isempty(ranks)
    error('No open contours found')
elseif length(ranks) < rank
    warning('Fewer than %i contours found. Setting rank = %i',rank,length(ranks))
    rank = length(ranks);
end
i = ranks(rank);
li = M(2,i);
cx = M(1,i+1:i+li);
cy = M(2,i+1:i+li);