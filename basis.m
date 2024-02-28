function b=basis(T,n,type)
% b=basis(T,n,type) makes basis functions of type 'hanning' or 'normal'. T
% is the length of the basis. n is how many to use. Typically T>>n, if the
% goal is to represent a filter with fewer parameters. Omit output argument
% to see how the basis functions look.

% (c) August 16, 2023 Lucas C Parra

r=T/n; 
switch type
    case 'hanning'
        for i=1:n b((1:round(r*4))+round((i-1)*r),i) = hanning(round(r*4)); end
        b=flipud(b(floor(r)+(1:T),:));
    case 'normal'
        t=(1:T)'; for i=1:n, b(:,i) = exp(-(t-(i-1)*r).^2/r^2); end
    otherwise
        b=eye(T);
end

if nargout==0, 
    plot([b sum(b,2)])
end


