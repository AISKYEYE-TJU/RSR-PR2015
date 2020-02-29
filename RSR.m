function W=RSR(X, lambda)
% regularized self-representation
%  min_X  || X - XW||_21 + lambda * ||W||_21

% Input : X n*m n number of samples m number of features
%         lambda  regularization parameter   
% Output: W m*m  self-representation matrix

iter_num = 50;
[n m] = size(X);

% initialize Gr and Gl
Gr = ones(m,1);
Gl = ones(n,1);

if m>n
    for iter = 1:iter_num
        G_R = spdiags(Gr,0,m,m);
        G_L = spdiags(Gl,0,n,n);
        
        % update W
        W = (G_R*X'*G_L*X+lambda*eye(m))\(G_R*X'*G_L*X);

        % update Gr
        wc = sqrt(sum(W.*W,2));
        Gr = 2*wc;
         
        % update Gl
        E = X*W-X;
        ec = sqrt(sum(E.*E,2)+eps);
        Gl = 0.5./ec;

        obj(iter) = sum(ec) + lambda*sum(wc);
    end;
else
    for iter = 1:iter_num
        G_R = spdiags(Gr,0,m,m);
        G_L = spdiags(Gl,0,n,n);
        
        % update W       
        W = G_R*X'*G_L*((X*G_R*X'*G_L+lambda*eye(n))\X);

        % update Gr
        wc = sqrt(sum(W.*W,2));
        Gr = 2*wc;

        % update Gl
        E = X*W-X;
        ec = sqrt(sum(E.*E,2)+eps);
        Gl = 0.5./ec;

        obj(iter) = sum(ec) + lambda*sum(wc);
    end;
end;