function [AA] = convertTo3d(A,sims)
    [M, N] = size(A);
    AA = zeros([sims.domain.Nx,sims.domain.Ny,sims.domain.Nz]);
    for n=1:N
        jj = n;
        kk = 0;
        for m=1:M
            if mod(m,sims.domain.Nx) == 0
                ii = sims.domain.Nx;
            elseif mod(m,sims.domain.Nx) == 1
                ii = 1;
                kk = kk + 1;
            else
                ii = mod(m,sims.domain.Nx);
            end
            AA(ii,jj,kk) = A(m,n);
        end
    end
end