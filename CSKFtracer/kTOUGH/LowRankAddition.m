    function [U,S,V]= LowRankAddition(U0,V0,U1,V1,r)
        % this function gives the low rank representation of A = U0V0' + U1V1'
        % so that A = U*S*V'
        % U and V
        U = [U0 U1]; %mxp, p= pu+pv
        V = [V0 V1]; %mxp
        % so that A + B = UV', with rank p, now compress to rank r
        %%%% Compute SVD of matrix product UV'%%%%%%%%%%%%%%%%%%%%
        % Compute QR decomposition of U and V, cost of (p^2(m+n))
        [Qu,Ru] = qr(U,0); % Ru has size of pxp
        [Qv,Rv] = qr(V,0); % Rv has size of pxp
        % SVD of RuRv' at O(p^3)
        [Ur,Sr,Vr] = svd(Ru*Rv');
        % SVD(U*V') = (Qu*Ur)*Sr*(Qv*Vr)'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        U = Qu*Ur(:,1:r);
        V = Qv*Vr(:,1:r);
        S = Sr(1:r,1:r);
    end