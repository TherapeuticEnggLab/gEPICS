%gEPICS: Generalized EPICS (Effective Pairwise Interactions for predicting Community Structure)

%The following code takes the number of species in the community,
%abundances in monoculture and leave-one-out subcommunities as inputs and
%produces effective pairwise interaction coeffiecients and the full-species abundance as
%outputs using gEPICS. Here, we assumed that the microbial community follows
%the generalized Lotka-Volterra model.
%We solve a system of linear equations (Ax=b) to obtain the effective pairwise interaction coefficients.
%The 'A'(Abundance_Matrix) and 'b'(The_Rate_matrix) are obtained as follows.
%We solve n^2 equations of which n*(n-1) equations correspond
%to abundances in leave-one-out subcommunities and n equations to
%monoculture abundances. The Abundance_Matrix is constructed in blocks. The
%first n blocks are of n-1 rows each, with each block corresponding to a
%leave-one-out experiment, starting from species 1 and ending with species n,
%a total of n*(n-1) equations. Within a block, each row corresponds to the abundance of each species present 
%in that leave-one-out experiment.The last block contains n equations each corresponding to 
%monoculture abundances. The effective pairwise interaction coefficients (unknowns), 
%are arranged in ascending order of species number(a11, a12, a13, ..., ann). 
%The right-hand side vector 'b' is derived from the intrinsic
%growth rates. We identified all the species that are extinct in leave-one-out
%subcommunities and eliminated those rows in Abundance_Matrix and The_Rate_matrix
%corresponding to these missing species, thereby making the system of equations
%under-determined. We use pseudo-inverse to solve the remaining equations and
%estimate all the unknowns.

clear all;
close all;
clc;
disp('Welcome to gEPICS');

%number of species in the community
num_species=input('\nEnter the number of species in the community=');


%monoculture abundances
mono_abun=input('\nEnter the monoculture abundances \n(Put your entries wrapped in square brackets. \nPut a space or a comma after each entry)=');

%leave-one-out abundances
loo_abun=input('\nEnter the abundances in leave-one-out subcommunities \n(Put your entries wrapped in square brackets. \n(Rows are subcommunities, columns are species)=');

global n %number of species in the community

n=num_species;

%intrinsic growth rate vector

r=ones(n,1);

Small_Matrix=zeros(n-1);
Large_Matrix=[]; 

%This command ensures that the diagonal elements of the leave-one-out abundances
%are zero
loo_abun(logical(eye(size(loo_abun)))) = 0;

%this loop creates matrix 'A'
for i=1:n
    %to create Medium_M, we first create Small_M 
    
    Small_Matrix=zeros(n-1,n);
    for j=1:n-1
        Small_Matrix(j,n*(j-1)+1:n*j)=loo_abun(i,:);
    end
    counter=1;
    flag=0;
    Medium_Matrix=zeros(n-1,n*n);
    
    %this loop combines Small_Matrix to create Medium_Matrix by adding zeros to
    %account for the missing species across leave-one-out cultures
    for k=1:n
        if k==i
            Medium_Matrix(:,(counter-1)*n+1:(counter)*n)=zeros(n-1,n);
            counter=counter+1;
            flag=1;
        
        else
            Medium_Matrix(:,(counter-1)*n+1:(counter)*n)=Small_Matrix(:,(counter-flag-1)*n+1:(counter-flag)*n);
            counter=counter+1;
        
        
    end
    end
    Large_Matrix=[Large_Matrix;Medium_Matrix]; 
end

Monoculture_appendment=zeros(n-1,n*n);

for i=1:n
    Monoculture_appendment(i,n*(i-1)+i)=mono_abun(1,i);
end

Abundance_Matrix=[Large_Matrix;Monoculture_appendment];

%this section of the code creates vector 'b'
Rate=[];
for i=1:n
    diagonal(i,1)=1;
    intrinsic_rate=r;
    intrinsic_rate(i)=[];
    Rate=[Rate;intrinsic_rate];
end

The_Rate_matrix=-1*[Rate;diagonal];

%this section of the code identifies the extinct species in any of the leave-one-out
%subcommunities. Abundances smaller than 10^-3 are negiligible and are
%considered to be extinct.
[subcommunity, species] = find(loo_abun<=10^-3);
Missing_Species = [subcommunity, species];
Missing_Species(subcommunity == species, :) = [];


%This section of the code makes the rows corresponding to the extinct species to zero in
%Larger_M matrix and The_big_R_matrix, thereby eliminating these equations
for i=1:height(Missing_Species)
Abundance_Matrix(((n-1)*(Missing_Species(i,1)-1))+Missing_Species(i,2),:)=zeros(1,n*n);
The_Rate_matrix(((n-1)*(Missing_Species(i,1)-1))+Missing_Species(i,2),:)=0;
end

%this command computes the effective pairwise interactions using pseudo-inverse 
intn_parameter=pinv(Abundance_Matrix)*The_Rate_matrix;

%the following commands rearrange effective pairwise interactions into matrix form
for i=1:n
    for j=1:n
        interaction_parameter(i,j)=intn_parameter(n*(i-1)+j);
    end
end

%this command predicts abundances in the full species community
full_species_abundance=-1*pinv(interaction_parameter)*r;

%Printing the outputs
disp('Effective Pairwise Interaction Coefficients estimated using gEPICS');
interaction_parameter
disp('Predicted Abundances in the full species community using gEPICS');
abundance = full_species_abundance'
