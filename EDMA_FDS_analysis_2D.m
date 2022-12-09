function EDMA_FDS_analysis_2D(res_path,nameOutput,id_cmp,age_range,num_samples,sample_step)
%% Initialization
if ~exist(res_path,'dir'), mkdir(res_path); end
clc; close all;

Label={'CNT';'Down';'Morquio';'NFM';'Noonan'};

addpath('utils_MATLAB/');
num_point   = 21;
scale       = 1;
idx         = [];
age         = [];
treatment   = [];
list_dir    = {};
sample_name = {};

% Schizophrenia
path_det        = 'Landmarks2D\Detection\';
path_dataset    = 'Landmarks2D\Dataset\';
l_dir       = dir(fullfile(path_dataset,'*.ini'));

for i=1:length(l_dir)
    list_dir = [list_dir; [path_det l_dir(i).name(1:end-3) 'txt']]; 
    fid  = fopen([path_dataset l_dir(i).name]);
	a = textscan(fid,'%s');
	fclose(fid);
    sample_name = [sample_name {l_dir(i).name(1:end-4)}];
    idx         = [idx str2num(a{1}{13})];
    age         = [age str2num(a{1}{7})];
    treatment   = [treatment str2num(a{1}{16})];
end 

id_age      = and(age>=min(age_range),age<max(age_range));
sample_name = sample_name(id_age);
list_dir    = list_dir(id_age);
age         = age(id_age);
idx         = idx(id_age);
treatment   = treatment(id_age);

id_id    = or(idx==0,idx==id_cmp);
sample_name = sample_name(id_id);
list_dir    = list_dir(id_id);
age         = age(id_id);
idx         = idx(id_id);
treatment   = treatment(id_id);

%% Read detections
x=[];
y=[];
puntos = [1:6,8:13,15,16,18:21];
for i=1:length(list_dir)
    [X, Y, ~] = read_det(list_dir{i}, num_point, false);
    x = [x; X(puntos)./scale];
    y = [y; Y(puntos)./scale];
end
num_point = length(puntos);

%% Calculate centroid
C = [mean(x,2), mean(y,2)];
CS = sum(((x-repmat(C(:,1),1,num_point)).^2+(y-repmat(C(:,2),1,num_point)).^2),2).^0.5;

%% Translation and scale centering
x = (x-repmat(C(:,1),1,num_point))./repmat(CS,1,num_point);
y = (y-repmat(C(:,2),1,num_point))./repmat(CS,1,num_point);

%% Calculate Form Matrix
[FM, id_ij] = FM_creator(x, y, num_point);

idx_CNT = idx==0;
idx_DIS = idx==id_cmp;

mean_dist_CNT = mean(FM(:,idx_CNT),2);
mean_dist_DIS = mean(FM(:,idx_DIS),2);

s_dist_CNT = var(FM(:,idx_CNT),0,2);
s_dist_DIS = var(FM(:,idx_DIS),0,2);
% 
n_CNT = sum(idx_CNT);
n_DIS = sum(idx_DIS);

diff_CNTvsDIS         = [];
diff_CNTvsDIS_All     = [];

%% Subgroups Analysis
n_Samples   = num_samples;
s_Samples   = sample_step;
n_SubGroups = 150;
IDX_SubGroups_CNT = [];
IDX_SubGroups_DIS = [];
CI = {};
for rep=1:n_SubGroups*(floor(n_Samples/s_Samples)+1)   
    n_DIS = (floor((rep-1)/n_SubGroups)*s_Samples);
    fprintf('Generating %i of %i Subgroups (%i Dis)\n',rep,n_SubGroups*(floor(n_Samples/s_Samples)+1),n_DIS);
    % Subgroups Generation
    [idx_SubGroupCNT,idx_SubGroupDIS] = Subgroups_Generation(n_Samples,n_DIS,idx_CNT,idx_DIS);
    IDX_SubGroups_CNT = [IDX_SubGroups_CNT;idx_SubGroupCNT];
    IDX_SubGroups_DIS = [IDX_SubGroups_DIS;idx_SubGroupDIS];
%     EDMA-II z test with montecarlo parametric bootstrapp
    alpha = 0.10;
    num_boot = 5000;
    CI_idx = floor([num_boot*(alpha/2), num_boot*(1-(alpha/2))]);

    mean_dist_SubGroup_CNT = mean(FM(:,idx_SubGroupCNT),2);
    mean_dist_SubGroup_DIS = mean(FM(:,idx_SubGroupDIS),2);
    cov_SubGroup_CNT = cov(FM(:,idx_SubGroupCNT)');
    cov_SubGroup_DIS = cov(FM(:,idx_SubGroupDIS)');
    cov_CNT=cov(FM(:,idx_CNT)');
    cov_DIS=cov(FM(:,idx_DIS)');

    mean_dist_CNT_boot=[];
    mean_dist_DIS_boot=[];
    mean_dist_SubGroup_CNT_boot=[];
    mean_dist_SubGroup_DIS_boot=[];
    
    for i = 1:num_boot
        mean_dist_CNT_boot  = [mean_dist_CNT_boot   ; mean(mvnrnd(mean_dist_CNT,cov_CNT,size(FM(:,idx_CNT),2)))];
        mean_dist_DIS_boot  = [mean_dist_DIS_boot   ; mean(mvnrnd(mean_dist_DIS,cov_DIS,size(FM(:,idx_DIS),2)))];
        mean_dist_SubGroup_CNT_boot  = [mean_dist_SubGroup_CNT_boot   ; mean(mvnrnd(mean_dist_SubGroup_CNT,cov_SubGroup_CNT,size(FM(:,idx_SubGroupCNT),2)))];
        mean_dist_SubGroup_DIS_boot  = [mean_dist_SubGroup_DIS_boot   ; mean(mvnrnd(mean_dist_SubGroup_DIS,cov_SubGroup_DIS,size(FM(:,idx_SubGroupDIS),2)))];
    end
    
    FDM_SubGroup_boot   = sort(mean_dist_SubGroup_CNT_boot - mean_dist_SubGroup_DIS_boot);
    FDM_All_boot        = sort(mean_dist_CNT_boot - mean_dist_DIS_boot);

    diff_CNTvsDIS       = [diff_CNTvsDIS; and(FDM_SubGroup_boot(CI_idx(1),:)<=0,FDM_SubGroup_boot(CI_idx(2),:)>=0)];
    diff_CNTvsDIS_All   = [diff_CNTvsDIS_All; and(FDM_All_boot(CI_idx(1),:)<=0,FDM_All_boot(CI_idx(2),:)>=0)];
    
    CI{end+1} = [FDM_SubGroup_boot(CI_idx(1),:)',FDM_SubGroup_boot(CI_idx(2),:)'];

end
p_All=sum(not(diff_CNTvsDIS_All),2)*100/size(id_ij,2);
p_SubGroup=sum(not(diff_CNTvsDIS),2)*100/size(id_ij,2);

figure(1)
for i=1:floor(n_Samples/s_Samples)+1
subplot(floor(n_Samples/s_Samples)+1,1,i)
histogram(p_SubGroup((n_SubGroups*(i-1))+1:n_SubGroups*i),0:2:100,'FaceColor',[0.2+(i*0.1) 0.2+(i*0.1) 0.2+(i*0.1)])
end
saveas(gcf,[res_path nameOutput],'fig')
saveas(gcf,[res_path nameOutput],'png')
end

function [x,y,z]=symetric_transform(x,y,z,Sag,Eyes,flip_id,Sym)
    for i=1:size(x,1)
        %plano sagital
        P=[x(i,:)' y(i,:)' z(i,:)'];
        P_Upper=mean(P(Sag.Upper,:),1);
        P_Lower=mean(P(Sag.Lower,:),1);
        Pm = repmat(mean([P_Upper;P_Lower],1),size(x,2),1);
        P=P-Pm;

        dx=P_Upper(1)-P_Lower(1);
        dy=P_Upper(2)-P_Lower(2);
        dz=P_Upper(3)-P_Lower(3);
        axy=atan2(dy,dx);
        azy=atan2(dy,dz);
        axz=atan2(dz,dx);
        Rx=[1 0 0; 0 cos(azy) sin(azy); 0 -sin(azy) cos(azy)];
        P=P*Rx;

        dx=P_Upper(1)-P_Lower(1);
        dy=P_Upper(2)-P_Lower(2);
        dz=P_Upper(3)-P_Lower(3);
        axy=atan2(dy,dx);
        azy=atan2(dy,dz);
        axz=atan2(dz,dx);
        Ry=[cos(axz) 0 -sin(axz); 0 1 0; sin(axz) 0 cos(axz)];
        P=P*Ry;

        dx=P_Upper(1)-P_Lower(1);
        dy=P_Upper(2)-P_Lower(2);
        dz=P_Upper(3)-P_Lower(3);
        axy=atan2(dy,dx);
        azy=atan2(dy,dz);
        axz=atan2(dz,dx);
        Rz=[cos(axy) -sin(axy) 0;sin(axy) cos(axy) 0; 0 0 1];
        P=P*Rz;
        x(i,:) = P(:,1)';
        y(i,:) = P(:,2)';
        z(i,:) = P(:,3)';

        %ojos
        P=[x(i,:)' y(i,:)' z(i,:)'];        
        P_EyeL=mean(P(Eyes.L,:),1);
        P_EyeR=mean(P(Eyes.R,:),1);
        Pm = repmat(mean([P_EyeL;P_EyeR],1),size(x,2),1);
        P=P-Pm;

        dx=P_EyeR(1)-P_EyeL(1);
        dy=P_EyeR(2)-P_EyeL(2);
        dz=P_EyeR(3)-P_EyeL(3);
        axy=atan2(dy,dx);
        azy=atan2(dy,dz);
        axz=atan2(dz,dx);
        Rx=[1 0 0; 0 cos(azy) sin(azy); 0 -sin(azy) cos(azy)];
        P=P*Rx;

        dx=P_EyeR(1)-P_EyeL(1);
        dy=P_EyeR(2)-P_EyeL(2);
        dz=P_EyeR(3)-P_EyeL(3);
        axy=atan2(dy,dx);
        azy=atan2(dy,dz);
        axz=atan2(dz,dx);
        Ry=[cos(axz) 0 -sin(axz); 0 1 0; sin(axz) 0 cos(axz)];
        P=P*Ry;

        dx=P_EyeR(1)-P_EyeL(1);
        dy=P_EyeR(2)-P_EyeL(2);
        dz=P_EyeR(3)-P_EyeL(3);
        axy=atan2(dy,dx);
        azy=atan2(dy,dz);
        axz=atan2(dz,dx);
        Rz=[cos(axy) -sin(axy) 0;sin(axy) cos(axy) 0; 0 0 1];
        P=P*Rz;
        x(i,:) = P(:,1)'-mean(P(:,1));
        y(i,:) = P(:,2)'-mean(P(:,2));
        z(i,:) = P(:,3)'-mean(P(:,3));
    end

    x_flip = -x(:,flip_id);
    y_flip = y(:,flip_id);
    z_flip = z(:,flip_id);

    x_sym  = (x+x_flip)./2;
    y_sym  = (y+y_flip)./2;
    z_sym  = (z+z_flip)./2;
    if Sym==1
        x = x_sym;
        y = y_sym;
        z = z_sym;
    end
end

function [FM, id_ij] = FM_creator(x, y, num_point)
    FM = zeros(sum(1:num_point-1),size(x,1));
    id_ij = [];
    for s=1:size(x,1)
        k = 1;
        for i=1:size(x,2)
            S_i = [x(s,i), y(s,i)];
            for j = i+1:size(x,2)
                S_j     = [x(s,j), y(s,j)];
                FM(k,s) = (sum((S_i-S_j).^2)).^0.5;
                if (s==1); id_ij   = [id_ij, [i;j]]; end;
                k       = k+1;
            end
        end
        FM(:,s) = FM(:,s)./geomean(FM(:,s));
    end      
end

function [idx_subCNT,idx_subDIS] = Subgroups_Generation(n_CNT,n_DIS,idx_CNT,idx_DIS)
CNT = find(idx_CNT);
DIS = find(idx_DIS);

new_CNT1 = datasample(CNT,n_CNT,'Replace',false);

idx_CNT(new_CNT1)=false;
CNT = find(idx_CNT);
new_CNT2 = datasample(CNT,n_CNT-n_DIS,'Replace',false);
new_DIS = datasample(DIS,n_DIS,'Replace',false);


idx_subCNT = false(size(idx_CNT));
idx_subCNT(new_CNT1) = true;
idx_subDIS = false(size(idx_DIS));
idx_subDIS([new_CNT2,new_DIS]) = true;
end

function [X, Y, Z] = read_det(name, num_points, DIM)
fdet=fopen(name,'r');
data=textscan(fdet,'%f');
fclose(fdet);
X   = [];
Y   = [];
Z   = [];
if (~DIM), data    = reshape(data{1},2,size(data{1},1)/2);
else data    = reshape(data{1},3,size(data{1},1)/3); end
X       = reshape(data(1,:),num_points,size(data,2)/num_points)';
Y       = reshape(data(2,:),num_points,size(data,2)/num_points)';
if (DIM), Z=reshape(data(3,:),num_points,size(data,2)/num_points)'; end
end

    