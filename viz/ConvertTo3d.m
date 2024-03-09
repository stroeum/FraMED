function [AA] = ConvertTo3d(A,B)
    [M, N] = size(A);
    AA = zeros(B');
    for n=1:N
        jj = n;
        kk = 0;
        for m=1:M
            if mod(m,B(1)) == 0
                ii = B(1);
            elseif mod(m,B(1)) == 1
                ii = 1;
                kk = kk + 1;
            else
                ii = mod(m,B(1));
            end
            AA(ii,jj,kk) = A(m,n);
        end
    end
end