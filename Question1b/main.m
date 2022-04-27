%% Solve equation -u''(x)=f(x) with the Dirichlet boundary condition 
clear all
close all
clc
%% Initial informations
ax=0.0;
bx=1.0;
cases=1;
N=16;% number of mesh points of first mesh
number_mesh=4;
number_mesh_point=zeros(number_mesh,1);
norm_max=zeros(number_mesh,1);
norm_l2=zeros(number_mesh,1);
norm_maxh1=zeros(number_mesh,1);
norm_h1=zeros(number_mesh,1);
%% Solve discrite solution and refine mesh

for inumber_mesh=1:number_mesh
    
    number_mesh_point(inumber_mesh)=N;
    
%% Create mesh point    
    x=zeros(N+1,1);
    for i=1:N+1
        x(i)=1-cos((i-1)*pi/(2*N));
    end
 %% Create matrix A   
    A=sparse(N-1,N-1);
%A = zeros(N-1,N-1);
    for i=1:N-1
        Delta1=x(i+2)-x(i+1);
        Delta2=x(i+1)-x(i);
        alpha=2/(Delta1*Delta2+Delta2^2);
        beta=Delta2/Delta1;
        gamma=Delta2/Delta1 + 1;
        if (i==1)
            A(i,i)=gamma;
            A(i,i+1)=-beta;
        elseif(i==N-1)
            A(i,i-1)=-1;
            A(i,i)=gamma;
        else
            A(i,i-1)=-1;
            A(i,i+1)=-beta;
            A(i,i)=gamma;
        end
        A(i,:)=alpha*A(i,:);
    end
%% Create vector b    
    b=zeros(N-1,1);
    for i=1:N-1
        b(i)=functionf(x(i+1),cases);
    end
%% Solve discrete solution
    u=A\b;
%% Get exact solution    
    u_ex=zeros(N+1,1);
    for i=1:N+1
        u_ex(i)=exact_solution(x(i),cases);
    end
%% Create discrete solution with boundary 
    u_dis=zeros(N+1,1);
    for i=1:N+1
        if (i==1)
            u_dis(i)=0;
        elseif(i==N+1)
            u_dis(i)=0;
        else
            u_dis(i)=u(i-1,1);
        end
    end
    
%% Calculate the error on L^infinity
    norm_max(inumber_mesh)=0.0;
    for i=1:N+1
        if (abs(u_dis(i)-u_ex(i)) > norm_max(inumber_mesh))
            norm_max(inumber_mesh)=abs(u_dis(i)-u_ex(i));
        end
    end
    
    norm_max(inumber_mesh)
%%  Calculate the error on L^2 

    norm_l2(inumber_mesh)=0;
    for i=1:N
        norm_l2(inumber_mesh)=norm_l2(inumber_mesh)+(u_dis(i)-u_ex(i))^2*(x(i+1)-x(i));
    end
    
    norm_l2(inumber_mesh)=(norm_l2(inumber_mesh))^(1/2);
    norm_l2(inumber_mesh)
%% Calculate the error on maxH1    

    norm_maxh1(inumber_mesh)=0;
    for i=1:N
        if (abs(((u_dis(i+1)-u_ex(i+1))-(u_dis(i)-u_ex(i)))/(x(i+1)-x(i))) > norm_maxh1(inumber_mesh))
            norm_maxh1(inumber_mesh)=abs(((u_dis(i+1)-u_ex(i+1))-(u_dis(i)-u_ex(i)))/(x(i+1)-x(i)));
        end
    end
    norm_maxh1(inumber_mesh)

%% Calculate the error on H1

    norm_h1(inumber_mesh)=0;
    for i=1:N
        norm_h1(inumber_mesh)=norm_h1(inumber_mesh)+(((u_dis(i+1)-u_ex(i+1))-(u_dis(i)-u_ex(i)))/(x(i+1)-x(i)))^2*(x(i+1)-x(i));
    end
    norm_h1(inumber_mesh)=(norm_h1(inumber_mesh))^(1/2);
%% Figure exact and dicrete solutions    
    figure
    plot(x,u_ex,'cyan h', x,u_dis,'magenta');
    xlabel('x');ylabel('value');
    title('Comparison between exact and discrete solutions with N=');
%     ylim([-0.02 0.12]);
    legend('Exact solution','Discrete solution');
    
%% Refine mesh (increse mesh point)    
    N=2*N;
end



%% Figure for errors respect to number of mesh point
figure
plot(log(number_mesh_point), -log(norm_max),'blue', log(number_mesh_point), -log(norm_l2), 'red',...
    log(number_mesh_point), -log(norm_maxh1), 'cyan', log(number_mesh_point), -log(norm_h1), 'magenta',...
    log(number_mesh_point), 2*log(number_mesh_point),'black');
xlabel('Log(MeshPoint)');ylabel('-Log(Error)');
title('Errors');
legend('norm_max','norm_l2','norm_maxh1','norm_h1','2x','Location','NorthEastOutside');  
