
% Copyright (C) 2024-2025 Sergii V. Siryk
% Copyright (C) 2024-2025 Vincenzo Di Florio
% 
% This program is free software: you can redistribute it and/or modify  
% it under the terms of the GNU General Public License as published by  
% the Free Software Foundation, either version 3 of the License, or  
% (at your option) any later version.  
% 
% This program is distributed in the hope that it will be useful,  
% but WITHOUT ANY WARRANTY; without even the implied warranty of  
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the  
% GNU General Public License for more details.  
% 
% You should have received a copy of the GNU General Public License  
% along with this program. If not, see <https://www.gnu.org/licenses/>. 

function [particles_coefficients,energy,linear_system]=multi_spheres_cg5(particles_params,medium_params,n_max,array_clebschgordan,cg_n_max,tol_linsolver, t_eps0_inv)
if numel(particles_params)<1
    error('At least one particle must be present. BYE.');
end
if nargin<6
    tol_linsolver=1e-10;
end
if test_overlapping(particles_params)==true
    error('SPHERES OVERLAP. BYE.');
end
if n_max>cg_n_max
    error('Must be n_max<=cg_n_max. BYE.');
end
kappa=medium_params.kappa;
epsilon_m=medium_params.epsilon;
Nnm=(n_max+1)^2;
linear_system.A=zeros(numel(particles_params)*Nnm);
linear_system.b=zeros([numel(particles_params)*Nnm,1]);

nmLM_mapping=zeros([n_max+1,2*n_max+1]);
inc=1;
for n=0:n_max
    nmLM_mapping(n+1,((-n):n)+(n_max+1))=inc:(inc+2*n);
    inc=inc+2*n+1;
end

if numel(particles_params)>1
    [XI_ij_matrices,Rij_mapping,sign_pattern_for_XI_Rji]=generate_XI_matrices(particles_params,n_max,kappa,nmLM_mapping,array_clebschgordan,cg_n_max);
end

%form left-hand side matrix A:
for i=1:numel(particles_params)
    tai=kappa*particles_params(i).radius;
    epsilon_i=particles_params(i).dielectric_constant;    
    linear_system.A(((i-1)*Nnm+1):i*Nnm,((i-1)*Nnm+1):i*Nnm)=full(generate_diagonal_matrix_n((0:n_max).*k_n(0:n_max,tai)*(epsilon_i-epsilon_m)/tai+epsilon_m*k_n((0:n_max)+1,tai),n_max));
    in_diag_i=generate_diagonal_matrix_n((0:n_max).*i_n(0:n_max,tai)*(epsilon_i-epsilon_m)/tai-epsilon_m*i_n((0:n_max)+1,tai),n_max);
    for j=1:(i-1)
        linear_system.A(((i-1)*Nnm+1):i*Nnm,((j-1)*Nnm+1):j*Nnm)=in_diag_i*(sign_pattern_for_XI_Rji.*XI_ij_matrices{Rij_mapping(j,i)});
    end
    for j=(i+1):numel(particles_params)
        linear_system.A(((i-1)*Nnm+1):i*Nnm,((j-1)*Nnm+1):j*Nnm)=in_diag_i*XI_ij_matrices{Rij_mapping(i,j)};
    end
end

%form right-hand side vector b:
for i=1:numel(particles_params)
    tai=kappa*particles_params(i).radius;
    epsilon_i=particles_params(i).dielectric_constant;
    vector_add_local=zeros([Nnm,1]);
    Hat_L00_i=sqrt(4*pi)*kappa*particles_params(i).charge/epsilon_i;
    n=0;
    vector_add_local(1)=(2*n+1)*epsilon_i*Hat_L00_i/(tai^(n+2));
    linear_system.b(((i-1)*Nnm+1):i*Nnm)=vector_add_local;
end

% [L,U]=ilu(sparse(linear_system.A),struct('type','ilutp','droptol',1e-1));
% coeffs=gmres(linear_system.A,linear_system.b,10,1e-10,100,L,U);
% clear L U;
Upsilon=zeros([size(linear_system.A,2),1]);
inc=numel(Upsilon);
for i=numel(particles_params):-1:1    
    Upsilon((inc-Nnm+1):inc)=generate_vec_temp(1./((2*(0:n_max)+1).*k_n(0:n_max,kappa*particles_params(i).radius)*particles_params(i).radius),n_max);
    inc=inc-Nnm;
end
Upsilon=spdiags(Upsilon,0,numel(Upsilon),numel(Upsilon));

linear_system.A=linear_system.A*Upsilon;

M=spdiags(diag(linear_system.A),0,size(linear_system.A,1),size(linear_system.A,2));
coeffs=gmres(linear_system.A,linear_system.b,10,tol_linsolver,100,M,[],M\linear_system.b);
coeffs=Upsilon*coeffs;
clear M Upsilon;

% coeffs=linear_system.A\linear_system.b;
% coeffs=full(sparse(linear_system.A)\sparse(linear_system.b));

% linear_system.A=A;
% linear_system.b=b;
% clear A b;

% inc=1;
% for i=1:numel(particles_params)
%     particles_coefficients(i).G_nm=coeffs(inc:(inc+Nnm-1));
%     inc=inc+Nnm;
% end
inc=numel(coeffs);
for i=numel(particles_params):-1:1
    particles_coefficients(i).G_nm=coeffs((inc-Nnm+1):inc);
    inc=inc-Nnm;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:numel(particles_params)
    tai=kappa*particles_params(i).radius;
    epsilon_i=particles_params(i).dielectric_constant;
    particles_coefficients(i).L_nm=zeros([Nnm,1]);
    vector_add_local=zeros([Nnm,1]);
    Hat_L00_i=sqrt(4*pi)*kappa*particles_params(i).charge/epsilon_i;
    n=0;
    vector_add_local(1)=-Hat_L00_i/(tai^(n+1));    
    temp=zeros([Nnm,1]);
    for j=1:(i-1)
        temp=temp+(sign_pattern_for_XI_Rji.*XI_ij_matrices{Rij_mapping(j,i)})*particles_coefficients(j).G_nm;
    end
    for j=(i+1):numel(particles_params)
        temp=temp+XI_ij_matrices{Rij_mapping(i,j)}*particles_coefficients(j).G_nm;
    end
    temp=generate_diagonal_matrix_n(i_n(0:n_max,tai),n_max)*temp;
    particles_coefficients(i).L_nm=generate_diagonal_matrix_n(k_n(0:n_max,tai),n_max)*particles_coefficients(i).G_nm+vector_add_local+temp;
    particles_coefficients(i).L_nm=generate_diagonal_matrix_n(1./(tai.^(0:n_max)),n_max)*particles_coefficients(i).L_nm;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
energy=0;
for i=1:numel(particles_params)
    energy=energy+particles_params(i).charge*particles_coefficients(i).L_nm(1)/sqrt(4*pi);
end
energy=t_eps0_inv*energy/2;
end

function res=k_n(n,x)
res=sqrt(2/pi)*besselk(n+1/2,x)./sqrt(x);
end

function res=i_n(n,x)
res=sqrt(pi/2)*besseli(n+1/2,x)./sqrt(x);
end

function res=test_overlapping(particles_params)
res=false;
for i=1:numel(particles_params)
    for j=(i+1):numel(particles_params)
        if norm(particles_params(j).center-particles_params(i).center)<=(particles_params(i).radius+particles_params(j).radius)
            res=true;
            break;
        end
    end
end
end

function [XI_ij_matrices,Rij_mapping,sign_pattern_for_XI_Rji]=generate_XI_matrices(particles_params,n_max,kappa,nmLM_mapping,array_clebschgordan,cg_n_max)
N_particles=numel(particles_params);
Rij_mapping=zeros(N_particles);
x=zeros([N_particles*(N_particles-1)/2,1]);
y=zeros(size(x));
z=zeros(size(x));
inc=1;
for i=1:N_particles
    for j=(i+1):N_particles
        x(inc)=particles_params(j).center(1)-particles_params(i).center(1);
        y(inc)=particles_params(j).center(2)-particles_params(i).center(2);
        z(inc)=particles_params(j).center(3)-particles_params(i).center(3);
        Rij_mapping(i,j)=inc;
        inc=inc+1;
    end
end
[phi,theta,Rij]=cart2sph(x(:),y(:),z(:));
theta=pi/2-theta;

knYnm=zeros([N_particles*(N_particles-1)/2,(2*n_max+1)*(2*n_max+2)/2]);
Ynm_mapping=zeros(2*n_max+1);
inc=1;
temp_exp=exp(1i*(0:(2*n_max)).*phi);
for n=0:(2*n_max)
    knYnm(:,inc:(inc+n))=k_n(n,kappa*Rij).*((((-1).^(0:n)).*(legendre(n,cos(theta),'norm')')/sqrt(2*pi)).*temp_exp(:,(0:n)+1));
    Ynm_mapping(n+1,(0:n)+1)=inc:(inc+n);
    inc=inc+n+1;
end
clear temp_exp;

% load array_cg10.mat array_clebschgordan cg_n_max;
offs_m1=cg_n_max+1;
offs_m2=2*cg_n_max+1;
offs_M=cg_n_max+1;

Nnm=(n_max+1)^2;
XI_ij_matrices=zeros([N_particles*(N_particles-1)/2,Nnm^2]);%zeros([(n_max+1)^2,(n_max+1)^2,N_particles*(N_particles-1)/2]);
sign_pattern_for_XI_Rji=ones((n_max+1)^2);
for n=0:n_max
    for m=(-n):0
        for L=0:n_max
            for M=(-L):L
                temp=0;
                m2=M-m;
                for l2=abs(n-L):(n+L)
                    if (abs(m2)<=l2)&&(mod(L+n+l2,2)==0)
                        if m2>=0
                            temp=temp+array_clebschgordan(n+1,0+offs_m1,l2+1,0+offs_m2,L+1,0+offs_M).*array_clebschgordan(n+1,m+offs_m1,l2+1,m2+offs_m2,L+1,M+offs_M).*sqrt(4*pi*(2*n+1).*(2*l2+1)./(2*L+1)).*knYnm(:,Ynm_mapping(l2+1,m2+1));
                        else
                            temp=temp+array_clebschgordan(n+1,0+offs_m1,l2+1,0+offs_m2,L+1,0+offs_M).*array_clebschgordan(n+1,m+offs_m1,l2+1,m2+offs_m2,L+1,M+offs_M).*sqrt(4*pi*(2*n+1).*(2*l2+1)./(2*L+1)).*conj(knYnm(:,Ynm_mapping(l2+1,-m2+1)))*((-1)^(-m2));
                        end
                    end
                end
                temp=((-1)^L)*temp;
                XI_ij_matrices(:,sub2ind([Nnm Nnm],nmLM_mapping(n+1,m+(n_max+1)),nmLM_mapping(L+1,M+(n_max+1))))=temp;
                XI_ij_matrices(:,sub2ind([Nnm Nnm],nmLM_mapping(n+1,-m+(n_max+1)),nmLM_mapping(L+1,-M+(n_max+1))))=((-1)^(m+M))*conj(temp);
                if mod(n+L,2)~=0
                    sign_pattern_for_XI_Rji(nmLM_mapping(n+1,m+(n_max+1)),nmLM_mapping(L+1,M+(n_max+1)))=-1;
                    sign_pattern_for_XI_Rji(nmLM_mapping(n+1,-m+(n_max+1)),nmLM_mapping(L+1,-M+(n_max+1)))=-1;
                end
            end
        end
    end
end
XI_ij_matrices=squeeze(num2cell(permute(reshape(XI_ij_matrices,[N_particles*(N_particles-1)/2,Nnm,Nnm]),[2 3 1]),[1 2]));
end

function res=generate_vec_temp(elements_on_diaglonal,n_max)
Nnm=(n_max+1)^2;
res=zeros([Nnm,1]);
inc=1;
for n=0:n_max
    res(inc:(inc+2*n))=elements_on_diaglonal(n+1);
    inc=inc+2*n+1;
end
end

function n_diag_matrix=generate_diagonal_matrix_n(elements_on_diaglonal,n_max)
Nnm=(n_max+1)^2;
n_diag_matrix=spdiags(generate_vec_temp(elements_on_diaglonal,n_max),0,Nnm,Nnm); %n_diag_matrix=sparse(diag(temp));
end