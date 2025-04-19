
clear
clc
tic
warning off

casename = 'WB2';
mips_opf_define

Y=[];
FOPFP=[];

F_num_file=0;
F_num_max=0;

NPS=50000;
for i=1:NPS
    theta = zeros(nb, 1);
    a = uu(nb+1: end) - ll(nb+1: end);
    b = ll(nb+1: end);
    VPQ = rand(nx-nb, 1);
    VPQ = a.* VPQ + b;   
    y0(:,i) = [theta; VPQ];
    x0=y0(:,i);
    [x, f, success, Output] = ...
        mips(f_fcn, x0, A, l, u, xmin, xmax, gh_fcn, hess_fcn,opt);
    if success==-1
        F_num_file=F_num_file+1;
    elseif success==0
        F_num_max=F_num_max+1;
    elseif success==1
        FOPFP=[FOPFP,x];
    end
end

plot(FOPFP(nb+1,:),FOPFP(nb+2,:),'r.')
axis([0.94 0.96 0.95 1.05])
t=toc;
fprintf('END of Point nummber is %d.Time is %8.3fs.',NPS,t)