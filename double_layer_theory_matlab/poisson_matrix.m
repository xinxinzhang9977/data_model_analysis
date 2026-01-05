function A = poisson_matrix(Nx, dx, eps, qs)

A = sparse(Nx,Nx);

for i = 2:Nx-1
    A(i,i-1) = eps/dx^2;
    A(i,i)   = -2*eps/dx^2;
    A(i,i+1) = eps/dx^2;
end

% x = 0 : Neumann BC (surface charge)
A(1,1) = -1;
A(1,2) = 1;

% x = L : reference potential
A(Nx,Nx) = 1;

end
