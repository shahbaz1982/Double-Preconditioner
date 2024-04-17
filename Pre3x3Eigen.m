% REFINED PRECONDITIONING STRATEGIES FOR DOUBLE
%SADDLE POINT PROBLEMS IN LIQUID CRYSTAL DIRECTOR
%MODELING by SHAHBAZ AHMAD

clear all
close all
clc
%Off-State problem is running. For On-state problem chage
% alpha =1.5 in Line 292.
global N
N=2^8; % size of the Matrix
K=9e-12;
[A,BB,D,CC]=CoMat(K,N);

[n nnA]=size(A);
[mC mmC]=size(CC);
[p ppD]=size(D);
B=BB';
C=CC';
[m mmB]=size(B);

Au = [A, B', C'; 
      B, zeros(m, m), zeros(m, p); 
      C, zeros(p, m), -D];

  % Construct the right-hand side vector b
f = rand(n, 1); % Example vector
g = rand(m, 1); % Example vector
h = rand(p, 1); % Example vector
b = [f; g; h];

% Define your system matrix A, preconditioners P and K, right-hand side vector b, and initial guess u0
SA=B*inv(A)*B';
SB=B*inv(A)*C';
SC=D+C*inv(A)*C';

tol = 1e-8; % tolerance for FGMRES

% Initialize FGMRES variables
[Ar Ac]=size(Au);
u0 = zeros(Ar, 1);
max_iter = 100; % maximum number of iterations for FGMRES

%---------------------1-PR1 Preconditioner---------------------------------------
Q1  = [eye(n,n), zeros(n, m) zeros(n, p); 
      B*inv(A), eye(m, m), zeros(m, p); 
      C*inv(A), zeros(p, m), eye(p, p)];
  
Q2  = [A, zeros(n, m) zeros(n, p); 
      zeros(m, n), SA, zeros(m, p); 
      zeros(p, n), zeros(p, m), SC];

Q3  = [eye(n,n), inv(A)*B', inv(A)*C'; 
      zeros(m, n), -eye(m, m), zeros(m, p); 
      zeros(p, n), zeros(p, m), -eye(p, p)];
  
PR1  = Q1*Q2*Q3;
KPR1=inv(Q3)*Q1';
tic  
[u1, flag1, relres1, iter1, resvec1] = RPCG(Au, PR1, KPR1, b, u0, tol, max_iter);

semilogy(resvec1,'-*');
disp('-----------------------------------------------------')
fprintf('Preconditioner     Flag       RES                 ITER\n')
disp('-----------------------------------------------------')
fprintf('%d                   %d        %d           %d\n',1,flag1,relres1,iter1)
toc
%----------------------2-PR2 Preconditioner---------------------------------------
Q4  = [A, zeros(n, m) zeros(n, p); 
      zeros(m, n), SA, zeros(m, p); 
      zeros(p, n), zeros(p, m), D];

PR2  = Q1*Q4*Q3;

KPR2=inv(Q3)*Q1';
tic  
[u2, flag2, relres2, iter2, resvec2] = RPCG(Au, PR2, KPR2, b, u0, tol, max_iter);

hold on
semilogy(resvec2,'-*')
 hold off
fprintf('%d                   %d        %d           %d\n',2,flag2,relres2,iter2);
toc
%-----------------------3-PMR1 Preconditioner---------------------------------------
Q5  = [eye(n,n), zeros(n, m) zeros(n, p); 
      B*inv(A), eye(m, m), zeros(m, p); 
      C*inv(A), SB'*inv(SA), eye(p, p)];
  
Q6  = [eye(n,n), inv(A)*B', inv(A)*C'; 
      zeros(m, n), -eye(m, m), -inv(SA)*SB; 
      zeros(p, n), zeros(p, m), -eye(p, p)];
  
PMR1  = Q5*Q2*Q6;
KPMR1=inv(Q6)*Q5';
tic
[u3, flag3, relres3, iter3, resvec3] = RPCG(Au, PMR1, KPMR1, b, u0, tol, max_iter);

hold on
semilogy(resvec3,'-o')
 hold off
fprintf('%d                   %d        %d           %d\n',3,flag3,relres3,iter3)
toc
%--------------------------4-PMR2 Preconditioner---------------------------------------

PMR2  = Q5*Q4*Q6;
KPMR2=inv(Q6)*Q5';
tic
[u4, flag4, relres4, iter4, resvec4] = RPCG(Au, PMR2, KPMR2, b, u0, tol, max_iter);

hold on
semilogy(resvec4,'-o')
 hold off
fprintf('%d                   %d        %d           %d\n',4,flag4,relres4,iter4)
toc

%---------------------------5-PAR1---------------------------------------------------
gam=100;
PAR1  = [A, 2*B', 2*C'; 
      B, SA, SB; 
      C, SB', -D+C*inv(A)*C'];


 PAR1N  = [eye(n,n), 2*inv(A)*B', 2*inv(A)*C'; 
      zeros(m, n), eye(m, m), zeros(m, p); 
      zeros(p, n), zeros(p, m), eye(p, p)];
 PAR1M  = [eye(n,n), zeros(n, m) zeros(n, p); 
      B*inv(A), eye(m, m), zeros(m, p); 
      C*inv(A), gam*eye(p, m), eye(p, p)];
  
KPAR1=inv(PAR1N)*PAR1M';
tic  
[u5, flag5, relres5, iter5, resvec5] = RPCG(Au, PAR1, KPAR1, b, u0, tol, max_iter);

hold on
semilogy(resvec5,'-+')
hold off
fprintf('%d                   %d        %d           %d\n',5,flag5,relres5,iter5)
toc
%---------------------------6-PAR2 Preconditioner------------------------------------------

PAR2  = [A, zeros(n, m) zeros(n, p); 
      B, -SA, -SB; 
      C, -SB', -SC];

PAR2N  = [eye(n,n), zeros(n, m) zeros(n, p); 
      zeros(m, n), eye(m, m), zeros(m, p); 
      zeros(p, n), zeros(p, m), eye(p, p)];
  
KPAR2=inv(PAR2N)*PAR1M';
tic
[u6, flag6, relres6, iter6, resvec6] = RPCG(Au, PAR2, KPAR2, b, u0, tol, max_iter);
fprintf('%d                   %d        %d           %d\n',6,flag6,relres6,iter6)
toc

hold on
semilogy(resvec6,'-+')
 
legend('PR1','PR2','PMR1','PMR2','PAR1','PAR2')
%title('Relative Residual Norms')
hold off

disp('-----------------------------------------------------')
disp('1-PR1,2-PR2,3-PMR1,4-PMR2,5-PAR1,6-PAR2')
disp('-----------------------------------------------------')
disp('Generating Residuals Graph')
disp ('Generating Eigenvalues Graphs ...')
%--------------------------------Eigenvalues Display------------------------------------------
figure;%A
eigenvalues_A = eig(full(Au));
plot((real(eigenvalues_A)),imag(eigenvalues_A),'o')
title('Eigenvalues of Matrix A ');
xlabel('Real Part');
ylabel('Imaginary Part');

figure;%PR1^{-1}A
eigenvalues_A = eig(full(inv(PR1)*Au));
plot((real(eigenvalues_A)),imag(eigenvalues_A),'o')
title('Eigenvalues of Matrix P_{R1}^{-1} A ');
xlabel('Real Part');
ylabel('Imaginary Part');

figure;%PR2^{-1}A
eigenvalues_A = eig(full(inv(PR2)*Au));
plot((real(eigenvalues_A)),imag(eigenvalues_A),'o')
title('Eigenvalues of Matrix P_{R2}^{-1} A ');
xlabel('Real Part');
ylabel('Imaginary Part');


figure;%PMR1^{-1}A
eigenvalues_A = eig(full(inv(PMR1)*Au));
plot((real(eigenvalues_A)),imag(eigenvalues_A),'o')
title('Eigenvalues of Matrix P_{MR1}^{-1} A ');
xlabel('Real Part');
ylabel('Imaginary Part');

figure;%PMR2^{-1}A
eigenvalues_A = eig(full(inv(PMR2)*Au));
plot((real(eigenvalues_A)),imag(eigenvalues_A),'o')
title('Eigenvalues of Matrix P_{MR2}^{-1} A ');
xlabel('Real Part');
ylabel('Imaginary Part');

figure;%PAR1^{-1}A
eigenvalues_A = eig(full(inv(PAR1)*Au));
plot((real(eigenvalues_A)),imag(eigenvalues_A),'o')
title('Eigenvalues of Matrix P_{AR1}^{-1} A ');
xlabel('Real Part');
ylabel('Imaginary Part');

figure;%PAR2^{-1}A
eigenvalues_A = eig(full(inv(PAR2)*Au));
plot((real(eigenvalues_A)),imag(eigenvalues_A),'o')
title('Eigenvalues of Matrix P_{AR2}^{-1} A ');
xlabel('Real Part');
ylabel('Imaginary Part');


%---------------------------------RPCG---------------------------------------------
function [u, flag1, relres1, iter1, resvec1] = RPCG(A, P, K, b, u0, tol, max_iter)
    % Input:
    % A: System matrix
    % P: Preconditioner for Pz = r
    % K: Preconditioner for Kv = z
    % b: Right-hand side vector
    % u0: Initial guess
    % tol: Tolerance for convergence
    % max_iter: Maximum number of iterations
    
    n = size(A, 1);
    m = size(P, 1);
    p = size(K, 1);
    
    u = u0;
    r = b - A * u;
    z = P \ r;
    p = z;
    v = K \ z;
    q = v;
    
    resvec1 = zeros(max_iter, 1);
    
    for iter = 1:max_iter
        alpha = -(v' * r) / (q' * A * p);
        u = u + alpha * p;
        r = r + alpha * A * p;
        
        z = P \ r;
        beta = -(v' * r) / (v' * K * q);
        p = z + beta * p;
        
        v = K \ z;
        q = v + beta * q;
        
        resvec1(iter) = norm(r);
        
        if norm(r) < tol
            break;
        end
    end
    
    % Check convergence
    if norm(r) <= tol
        flag1 = 0; % Convergence achieved
    else
        flag1 = 1; % Maximum number of iterations reached without convergence
    end
    
    relres1 = norm(r) / norm(b); % Relative residual
    iter1 = iter; % Number of iterations
end


%--------------------------Function----------------------------
function [A,B,C,D]=CoMat(K,N)
%clear all
%format short e

global n Np1 dz dzi alsq beta 
global u v w phi lambda 
global precon solver msolve

% LC parameters
%K=9e-12;
eps_0=8.8542e-12; 
eps_par=15; eps_perp=5; eps_a=eps_par-eps_perp;
alphac=sqrt(3)*pi/2; % 2.721
beta=eps_perp/eps_a;

% set desired alpha
alpha=0.5*alphac; % no electric field for Off-state
%alpha=1.5*alphac; % On-state
alsq=alpha^2;

% find resulting voltage
V=sqrt((K*alsq)/(eps_0*eps_a));

fprintf('alpha: %8.4e, beta: %8.4e, voltage: %8.4e\n', ...
         alpha,beta,V)
 
% set grid size N
%N=2^2;
fprintf('\n N=%3i\n', N);

% construct coefficient matrix
solver=0; % 0 for full matrix, 1 for reduced matrix
[H,A,B,C,D]= mgen(V);
end

function [H,A,B,C,D]=mgen(V)
% constructs initial matrix for 1D TN test problem

global n N Np1 dz dzi alsq beta 
global u v w phi lambda
global precon solver msolve

% number of unknowns
n=N-1; Np1=N+1;

% discretisation parameter
dz=1/N; dzi=1/dz;

% coordinates 
x=ones(Np1,1); y=ones(Np1,1); z=linspace(0,1,Np1)';

% choose initial guess 
init=2;
if init==1 
   % pure twist
   angle=0; theta=angle*ones(Np1,1); psi=linspace(0,pi/2,Np1)';
   u=cos(theta).*cos(psi); v=cos(theta).*sin(psi); w=sin(theta);
   phi=linspace(0,1,Np1)'; 
elseif init==2
   % twist and tilt
   theta=sin(pi*z);  psi=linspace(0,pi/2,Np1)';  % top curve 
%  theta=sin(-pi*z); psi=linspace(0,-pi/2,Np1)'; % bottom curve 
   u=cos(theta).*cos(psi); v=cos(theta).*sin(psi); w=sin(theta);
   phi=linspace(0,1,Np1)'; 
end

% inital guess for lambda
initl=2;
if initl==1 
   % pure twist
   lambda=-pi^2/4*dz*ones(Np1,1); lambda(1)=0; lambda(end)=0;
elseif initl==2
   % twist and tilt
   lambda=dz*alsq*ones(Np1,1); lambda(1)=0; lambda(end)=0;
end 

% apply boundary conditions
u(1)=1;v(1)=0;w(1)=0;u(end)=0;v(end)=1;w(end)=0;  % top curve
%u(1)=1;v(1)=0;w(1)=0;u(end)=0;v(end)=-1;w(end)=0; % bottom curve
phi(1)=0; phi(end)=1; 

% construct full solution vector
nvec=[u; v; w];
sol=[nvec; lambda; phi];

% construct Hessian
[H,A,B,C,D]=set_hessian;

if solver==1
   % construct nullspace matrix
   val=sqrt(v(2:N).^2+w(2:N).^2);
   Nmat(1:n,n+1:2*n)=spdiags(val,0,n,n);
   Nmat(n+1:2*n,1:n)=spdiags(-w(2:N)./val,0,n,n);
   Nmat(n+1:2*n,n+1:2*n)=spdiags(-(v(2:N).*u(2:N))./val,0,n,n);
   Nmat(2*n+1:3*n,1:n)=spdiags(v(2:N)./val,0,n,n);
   Nmat(2*n+1:3*n,n+1:2*n)=spdiags(-(w(2:N).*u(2:N))./val,0,n,n);
   % construct rhs
   [grad,grad_n,grad_l,grad_p]=setup_gradient; rhs=-grad;
   % set up reduced system without B
   xhat=-B*((B'*B)\grad_l);
   temp1=Nmat'*D;
   temp2=Nmat'*A*Nmat;
   H=[temp2 temp1; temp1' -C];
end


end

function [Hhat,Ahat,Bhat,Chat,Dhat]=set_hessian
% sets up Hessian with boundary equations ignored

global n N Np1 dz dzi alsq beta 
global u v w phi lambda

os=ones(n,1); 

% set up midpoint values
wsqhalf=(w(2:Np1).^2+w(1:N).^2)/2;
ahalf=alsq*(beta+wsqhalf);
phidiff=phi(2:Np1)-phi(1:N);

% set up Ahat
vc=2+dz*lambda(2:N); 
vr=-os; vl=-os;
A11=dzi*spdiags([vl vc vr],-1:1,n,n);
A22=A11;
vp=phidiff(1:n).^2+phidiff(2:N).^2;
A33=A11-dzi*alsq/2*spdiags(vp,0,n,n);
%Ahat=blkdiag(A11,A22,A33);
Ahat(1:n,1:n)=A11;
Ahat(n+1:2*n,n+1:2*n)=A22;
Ahat(2*n+1:3*n,2*n+1:3*n)=A33;

% set up Bhat
BU=spdiags(u(2:N),0,n,n); 
BV=spdiags(v(2:N),0,n,n); 
BW=spdiags(w(2:N),0,n,n); 
%Bhat=[BU; BV; BW];
Bhat(1:n,1:n)=BU;
Bhat(n+1:2*n,1:n)=BV;
Bhat(2*n+1:3*n,1:n)=BW;

% set up Chat
vc=ahalf(1:n)+ahalf(2:N);
vl=[-ahalf(2:n); 0]; vr=[0; -ahalf(2:n)]; 
Chat=dzi*spdiags([vl vc vr],-1:1,n,n);

%keyboard
% set up Dhat
mul=alsq*dzi*spdiags(w(2:N),0,n,n);
vc=phi(1:n)-2*phi(2:N)+phi(3:Np1);
vl=[phidiff(2:n);0]; vr=[0;-phidiff(2:n)];
DW=mul*spdiags([vl vc vr],-1:1,n,n);
Dhat=spalloc(3*n,n,3*n);
Dhat(2*n+1:3*n,:)=DW;

%e=ones(n,1); z=zeros(n,1); mat=alsq*spdiags([e z -e],-1:1,n,n);

%ev=eig(full(inv(Chat)));
%nC=norm(full(inv(Chat)));
%disp([nC nC*dz nC/dz])

% set up Hhat
temp=[Bhat Dhat];
tC=spalloc(2*n,2*n,3*n);
tC(n+1:2*n,n+1:2*n)=-Chat;
%Hhat=[Ahat temp; temp' tC];
Hhat(1:3*n,1:3*n)=Ahat;
Hhat(1:3*n,3*n+1:5*n)=temp;
Hhat(3*n+1:5*n,1:3*n)=temp';
Hhat(3*n+1:5*n,3*n+1:5*n)=tC;


end
