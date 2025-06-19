clear all, close all

kmax = 

% Step 1 : generating the tessellation
l = [-1,1;-1,1];
Nv = 120;
l_A = .1*l;
tesstype = 'regular';
[C0,v0,dA0,~,Nv0,~,~] = TessellationGenerator(l,Nv,l_A,tesstype);
[v0,C0,Nv0,Nc0,dX0,Tess0,ext_ind0] = TessellationFinder(Nv0,v0,C0,l_A,'regular');

figure
triplot(C0, v0(1,:), v0(2,:), Color=[0.3010 0.7450 0.9330], LineWidth=0.2)
xlabel('x1')
ylabel('x2')
grid on;
figure
trisurf(C0,v0(1,:),v0(2,:),zeros(length(v0),1),'FaceAlpha',.6)
grid on;


% Step 2 : generating data
Nd = 300;
rt = 1.1;
system = 1; % 0,1,2,3
nrm = 2; % 1,2,Inf
M0 = LC(rt*l(1,2),func(system),nrm);
dtype = 1; % 0,1,2
[x0,f0,Nd0] = DataGenerator(Nd,rt*l,v0,Tess0,system,M0,nrm,dtype);

figure
triplot(C0, v0(1,:), v0(2,:), Color=[0.3010 0.7450 0.9330], LineWidth=0.2)
hold on
scatter(x0(1,:), x0(2,:), 'MarkerEdgeColor', [0.5 0.4470 0.7410])
quiver(x0(1,:), x0(2,:), f0(1,:), f0(2,:), Color=[0.8500 0.3250 0.0980])
grid on
xlabel('\theta')
ylabel('\omega')


% Step 3 : optimisation problem
eps = 1e-2; % negativity tolerance
eta = 5; % slack margin
alpha = 200; % upper limit of the Lyapunov function
beta = 0; % lower limit of the Lyapunov function
lim_mode = 0; % additional constraints on function edges (0:none,1:inner,2:outer,3:both)
clear gs0 bs0 vals0 ss0 sum_ss0
[gs0,bs0,vals0,ss0,sum_ss0] = OptimisationProblem(x0,f0,M0,v0,C0,eps,eta,alpha,beta,dA0,dX0,lim_mode,nrm,0,0,0);

figure
trisurf(C0,v0(1,:),v0(2,:),vals0,'FaceAlpha',0.6)
xlabel('\theta')
ylabel('\omega')
zlabel('V(x)')
grid on
colorbar

fprintf('Computed optimal value : %f\n', -eta*length(C0)*3)
fprintf('Obtained optimal value : %f\n', sum_ss0)
fprintf('Error : %f\n', sum_ss0 + eta*length(C0)*3)


% Step 4 : validating the results
clear bad_ind0 bad_triangles0 bad_points0
bad_ind0 = [];
for c = 1: length(ss0)
    for i = 1:3
        if ss0(i,c) >= 0
            bad_ind0 = [bad_ind0; [c, i, ss0(i,c)]];
        end
    end
end

figure
if not(size(bad_ind0)==0)
    bad_points0 = zeros(2, length(bad_ind0(:,1)));
    bad_triangles0 = zeros(length(bad_ind0(:,1)), 3);
    for ml = 1:length(bad_ind0(:,1))
        bad_points0(:,ml) = [v0(1,C0(bad_ind0(ml,1),bad_ind0(ml,2))); v0(2,C0(bad_ind0(ml,1),bad_ind0(ml,2)))]; % points' coordinates
        bad_triangles0(ml,:) = C0(bad_ind0(ml,1),:); % triangles' indices
    end
    scatter(bad_points0(1,:), bad_points0(2,:), 50, 'yellow', 'filled')
    hold on
    scatter(bad_points0(1,:), bad_points0(2,:), bad_ind0(:,3).*200./max(bad_ind0(:,3)), 'r', 'filled')
    triplot(bad_triangles0, v0(1,:), v0(2,:))
end
scatter(v0(1,:), v0(2,:), 3, 'b')
grid on

[pointlist0,a0] = sublevelset(dX0,C0,v0,gs0,bs0,vals0,1,[.2 .6 1]);

%% 
clear v1 C1 Tess1
[v1,C1,Nv1,Nc1,dX1,Tess1,ext_ind1] = TessellationFinderV2(250,v0,bad_triangles0,0.8*l_A,'random',pointlist0,0.025);

figure
triplot(C_new, v_new(1,:), v_new(2,:), Color=[0.3010 0.7450 0.9330], LineWidth=0.2)
xlabel('x1')
ylabel('x2')
grid on;
figure
trisurf(C_new,v_new(1,:),v_new(2,:),zeros(length(v_new),1),'FaceAlpha',.6)
grid on

l_new = [min(v_new(1,:)),max(v_new(1,:));min(v_new(2,:)),max(v_new(2,:))];
[x_new,f_new,Nd_new] = DataGenerator(350,rt,l_new,v_new,Tess_new,system,M,nrm,1);

figure
triplot(C_new, v_new(1,:), v_new(2,:), Color=[0.3010 0.7450 0.9330], LineWidth=0.2)
hold on
scatter(x_new(1,:), x_new(2,:), 'MarkerEdgeColor', [0.5 0.4470 0.7410])
quiver(x_new(1,:), x_new(2,:), f_new(1,:), f_new(2,:), Color=[0.8500 0.3250 0.0980])
grid on
xlabel('\theta')
ylabel('\omega')

clear gs_new bs_new vals_new ss_new sum_ss_new
[gs_new,bs_new,vals_new,ss_new,sum_ss_new] = OptimisationProblem(x_new,f_new,M,v_new,C_new,eps,eta,alpha,beta,dA,dX_new,lim_mode,nrm,0,0,0);
fprintf('Computed optimal value : %f\n', -eta*length(C_new)*3)
fprintf('Obtained optimal value : %f\n', sum_ss_new)
fprintf('Error : %f\n', sum_ss_new + eta*length(C_new)*3)

figure
trisurf(C_new,v_new(1,:),v_new(2,:),vals_new,'FaceAlpha',0.6)
xlabel('\theta')
ylabel('\omega')
zlabel('V(x)')
grid on
colorbar

clear bad_ind_new bad_triangles_new bad_points_new
bad_ind_new = [];
for c = 1: length(ss_new)
    for i = 1:3
        if ss_new(i,c) >= 0
            bad_ind_new = [bad_ind_new; [c, i, ss_new(i,c)]];
        end
    end
end