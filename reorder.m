function out = reorder(in)
% This function will get sequence of labels and rename them so the most
% repeted label be 1 then 2 and so on ...

nLabels = 5;
labels_repeat = zeros(1, nLabels);

for i = 1 : length(in)
    labels_repeat(in(i)) = labels_repeat(in(i)) + 1;
end

[labels_repeat_new, I] = sort(labels_repeat);
out = in;

for i = 1 : length(in)
    out(i) = find(I == in(i));
end

end

