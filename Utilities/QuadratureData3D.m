function [ NIps , Weights , Coords ] = QuadratureData3D( Type )

% NIps = Number of Ips
% Weights = Local Weights for each ip
% Coords  = Local Coordinates for each ip

switch Type
    case "Gauss 1x1x1"
        NIps = 1; 
        Weights = 8.0;
        Coords  = [ 0.0 , 0.0 , 0.0 ]; %Ip1
        
    case "Gauss 2x2x2"
        NIps = 8;
        Weights = 1.0 * ones(1,8);
        Coords  = (3^(-1/2)) * [-1, -1, -1;     %ip1
                                 1, -1, -1;     %ip2
                                 1,  1, -1;     %ip3
                                -1,  1, -1;     %ip4
                                -1, -1,  1;     %ip5
                                 1, -1,  1;     %ip6
                                 1,  1,  1;     %ip7
                                -1,  1,  1];    %ip8
    
    case "Gauss 3x3x3"
        NIps = 27;

        posgl = sqrt(0.6)*[-1 0 1];
        weigl = [5 8 5]/9;
        nlocs = 3;
        
        %Coords = [reshape(repmat(posgl',nlocs^2,1),[],1)  repmat([reshape(repmat(posgl',nlocs,1),[],1) reshape(repmat(posgl,nlocs,1),[],1)],nlocs,1)];
        Coords = [repmat([reshape(repmat(posgl',nlocs,1),[],1) reshape(repmat(posgl,nlocs,1),[],1)],nlocs,1)  reshape(repmat(posgl,nlocs^2,1),[],1)];
        Weights = reshape(reshape(weigl'*weigl,[],1)*weigl,1,[]);

    case "Gauss 5x5x5"
        NIps = 125;

        posgl = [ -sqrt(5+2*sqrt(10/7))/3 , -sqrt(5-2*sqrt(10/7))/3 , 0 , sqrt(5-2*sqrt(10/7))/3 , sqrt(5+2*sqrt(10/7))/3 ];
        weigl = [ (322-13*sqrt(70))/900 , (322+13*sqrt(70))/900 , 128/225 , (322+13*sqrt(70))/900 , (322-13*sqrt(70))/900];
        nlocs = 5;

        Coords = [repmat([reshape(repmat(posgl',nlocs,1),[],1) reshape(repmat(posgl,nlocs,1),[],1)],nlocs,1)  reshape(repmat(posgl,nlocs^2,1),[],1)];
        Weights = reshape(reshape(weigl'*weigl,[],1)*weigl,1,[]);

    otherwise 
        error('Quadrature option not available. Choose between "Gauss 1x1x1" and "Gauss 2x2x2"');        
end