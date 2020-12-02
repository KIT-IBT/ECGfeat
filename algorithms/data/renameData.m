clear
close all
d=dir('.')
d2=d([d.isdir]~=1);
s=string({d2.name});
ind1=contains(s,'.txt');
ind2=contains(s,'ecgout');
ind2change=find(ind1&ind2);
for j=1:1:length(ind2change)
    copyfile(s(ind2change(j)),erase(s(ind2change(j)),"ecgout_"))
end
