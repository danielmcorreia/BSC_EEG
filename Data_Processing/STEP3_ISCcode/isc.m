function [ISC,W, A, ISC_persubject, ISC_persubject_sum, ISC_perwindow] = isc(X, gamma, fs, Nsec, overlap)

    % T samples, D channels, N subjects
    [T,D,N] = size(X);

    % now start the ISC code

    % compute cross-covariance between all subjects i and j
    Rij = permute(reshape(cov(X(:,:)),[D N  D N]),[1 3 2 4]);

    % compute within- and between-subject covariances
    Rw = 1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
    Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects

    % shrinkage regularization of Rw
    Rw_reg = (1-gamma)*Rw + gamma*mean(eig(Rw))*eye(size(Rw));

    % compute correlated components W using regularized Rw, sort components by ISC
    [W,ISC]=eig(Rb,Rw_reg); [ISC,indx]=sort(diag(ISC),'descend'); W=W(:,indx);
    
    %1. [V,D] = eig(A) returns diagonal matrix D of eigenvalues and matrix V whose columns are the corresponding right eigenvectors
    %2. The component projections that capture the largest correlation
    %between subjects (ISC) are the eigenvectors vi of matrix 𝑹𝑤−1𝑹𝑏 with
    %the strongest eigenvalues
    %3. Spatial filter eigenvalues serve as component coefficients.
    %4. Conc.: ISC here are the eigenvalues of the spatial filters and it's
    %already a ISC value (it's similar to the average of the ISC per
    %subject, because the latest is calculated by comparing each subject's
    %activity with this "general" one)

    A=Rw*W/(W'*Rw*W);

    % Compute ISC resolved by subject, see Cohen et al.
    % Every row corresponds to one component, every column corresponds to each subject.

    % ISC persubject   

    for i=1:N
        Rw=0; for j=1:N, if i~=j, Rw = Rw+1/(N-1)*(Rij(:,:,i,i)+Rij(:,:,j,j)); end; end
        Rb=0; for j=1:N, if i~=j, Rb = Rb+1/(N-1)*(Rij(:,:,i,j)+Rij(:,:,j,i)); end; end
        ISC_persubject(:,i) = diag(W'*Rb*W)./diag(W'*Rw*W);
    end

    % Sum the 3 strongest components and save it as a vector
    ISC_persubject_sum=sum(ISC_persubject(1:3,:),1);

    % ISC across time windows 
    
    t = 1:floor((Nsec-overlap)*fs+1):(T-fs*Nsec);
    for tau=1:length(t)
        Xt = X(t(tau):(t(tau)+Nsec*fs-1),:,:);
        Rij = permute(reshape(cov(Xt(:,:)),[D N  D N]),[1 3 2 4]);
        Rw =  1/N* sum(Rij(:,:,1:N+1:N*N),3);  % pooled over all subjects
        Rb = 1/(N-1)/N*(sum(Rij(:,:,:),3) - N*Rw);  % pooled over all pairs of subjects
        ISC_perwindow(:,tau) = diag(W'*Rb*W)./diag(W'*Rw*W);
    end