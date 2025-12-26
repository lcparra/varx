function [pvalue,stat]=similarity_shuffle_test(M1,M2,Nrand)
% p=similarity_shuffle_test(M1,M2,Nrand) Shuffle test to determine if two
% matrices (or vectors) are similar. Nrand if the number of shuffles. The
% test shuffles elements in each matrix at random and then computes as test
% statistic (stat) the symmetrized Tanimoto distance between the original and
% shuffled matrices. A small p-value indicates that the the matrices are
% similar. The minimal p-value is set to 1/Nrand if no shuffle was larger
% than the original test statistic. If no putput argument is given, the
% matrices and null ditribution for the test statistic are shown.
% 
% The method was developed by Behtash Babadi and is explained in 
% https://doi.org/10.64898/2025.12.22.696055 

% (c) December 22, 2025 Lucas C Parra 
%                       based on code from Behtash Babadi and Jens Mdsen

if nargin<3, Nrand = 10000; end

if nargout<1
end

M1r = M1(:) / norm(M1(:),"fro");
M2r = M2(:) / norm(M2(:),"fro");

tanimoto = @(a,b) (a'*b) / ( (a'*a) + (b'*b) - (a'*b) );

for i=1:Nrand+1
    stat(i) = tanimoto(M2r,M1(:))/2 + tanimoto(M1r,M2(:))/2; % symmetrized Tanimoto distance
    M1r = M1r(randperm(numel(M1))); 
    M2r = M2r(randperm(numel(M2))); 
end
pvalue = max(mean(stat(2:end)>stat(1)),1/Nrand);

if nargout<1
    clims = [min([M1r;M2r]), max([M1r;M2r])];
    subplot(2,2,1); imagesc(M1); title('M_1'); axis square; clim(clims); colorbar;
    subplot(2,2,2); imagesc(M2); title('M_2'); axis square; clim(clims); colorbar;
    subplot(2,1,2);
    hist(stat(2:end));hold on; ax=axis;
    plot(stat(1)*[1 1],ax(3:4),'r'); hold off
    legend({'shuffle',['p=' num2str(pvalue,2)]})
    title(['Shuffle elements of between M_1 and M_2']);  
    xlabel('SymmetrizedTanimoto(M_1,M_2)')
else
    stat = stat(1);
end



return

% some test code

% some arbitrary correlation matrix to test
R=corr(rand(10));  
% look the same and the test says they are similar (p = 0)
similarity_shuffle_test(R+0.1*rand(10),R+0.1*rand(10))
% look very different and the test says they are not significantly similar (p = 1)
similarity_shuffle_test(R+10*rand(10),R+10*rand(10)) 
% look the same and the test says they are similar (p = 0)
similarity_shuffle_test(R,R+0.1) 