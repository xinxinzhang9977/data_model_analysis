function [Ws, Wv] = build_weights(x, R)

Nx = length(x);
dx = x(2)-x(1);

Ws = zeros(Nx);
Wv = zeros(Nx);

for i = 1:Nx
    for j = 1:Nx
        d = abs(x(i)-x(j));

        % --- surface (charge) weight ---
        if d <= R
            Ws(i,j) = 1/(2*R);
        end

        % --- volume (packing) weight ---
        if d <= R
            Wv(i,j) = 3*(R^2-d^2)/(4*R^3);
        end
    end
end

Ws = Ws * dx;
Wv = Wv * dx;

end
