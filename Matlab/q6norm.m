%function output=q6norm(data,radius,scalex,scaley,scalez);
function output = q6norm(directory,filestem,fileend,outstem,radius,number,scalex,scaley,scalez); 
	% start from particle position list, radius for neighbors, file number
    % outputs 4 column vector with position and # of Q6 neighbors of each particle
    
   fle = [directory filestem int2str(number) fileend];
   data = load(fle);
   positions = data(:,1:3);
   [numpart, b] = size(positions);
   positions(:,1) = positions(:,1)*scalex;
   positions(:,2) = positions(:,2)*scaley;
   positions(:,3) = positions(:,3)*scalez;

%    Q6allbonds = zeros(numpart,dumber);
    Q6bonds = zeros(numpart,4);
    XRqmIqm = zeros(numpart,29);
    X = positions(:,1);
    Y = positions(:,2);
    Z = positions(:,3);
    Q6bonds(:,1) = X;
    Q6bonds(:,2) = Y;
    Q6bonds(:,3) = Z;
    
    for i = 1:numpart
        xi = positions(i,1);
        yi = positions(i,2);
        zi = positions(i,3);
        R = sqrt((X-xi).*(X-xi) + (Y-yi).*(Y-yi) + (Z-zi).*(Z-zi));
        radii = [X Y Z R];
        Neighbor = radii;
        Neighbor(Neighbor(:,4)>radius,:)=[];
        Neighbor(Neighbor(:,4)==0,:)=[];
        [NN,a]=size(Neighbor);            
        RYi = zeros(NN,13);
        IYi = zeros(NN,13);
        
        for j = 1:NN
            xj = Neighbor(j,1);
            yj = Neighbor(j,2);
            zj = Neighbor(j,3);
            rij = Neighbor(j,4);
            
            costh = (zj-zi)/rij;
            sinth = sqrt(1 - costh*costh);
            if isequal(sinth,0)
                cosphi = 0;
                sinphi = 0;
            else 
                cosphi = (xj-xi)/(rij*sinth);
                sinphi = (yj-yi)/(rij*sinth);
            end
            phi = acos(cosphi);
            
            d0 = (1-costh^2)^6;
            d1 = -12*(1-costh^2)^5*costh;
            d2 = 120*(1-costh^2)^4*costh^2-12*(1-costh^2)^5;
            d3 = -960*(1-costh^2)^3*costh^3+360*(1-costh^2)^4*costh;
            d4 = 5760*(1-costh^2)^2*costh^4-5760*(1-costh^2)^3*costh^2+360*(1-costh^2)^4;
            d5 = -23040*(1-costh^2)*costh^5+57600*(1-costh^2)^2*costh^3-14400*(1-costh^2)^3*costh;
            d6 = 46080*costh^6-345600*(1-costh^2)*costh^4+259200*(1-costh^2)^2*costh^2-14400*(1-costh^2)^3;
            d7 = 967680*costh^5-2419200*(1-costh^2)*costh^3+604800*(1-costh^2)^2*costh;
            d8 = 9676800*costh^4-9676800*(1-costh^2)*costh^2+604800*(1-costh^2)^2;
            d9 = 58060800*costh^3-21772800*(1-costh^2)*costh;
            d10 = 239500800*costh^2-21772800;
            d11 = 479001600*costh;
            d12 = 479001600;
             
            for m=-6:6
                RYi(j,m+7) = 0.00002207264*((-1)^(6+m))*sqrt(myfactorial(6-m)/myfactorial(6+m))*cos(m*phi)*sinth^m;
                IYi(j,m+7) = 0.00002207264*((-1)^(6+m))*sqrt(myfactorial(6-m)/myfactorial(6+m))*sin(m*phi)*sinth^m;
            end
            RYi(j,1) = RYi(j,1)*d0;
            RYi(j,2) = RYi(j,2)*d1;
            RYi(j,3) = RYi(j,3)*d2;
            RYi(j,4) = RYi(j,4)*d3;
            RYi(j,5) = RYi(j,5)*d4;
            RYi(j,6) = RYi(j,6)*d5;
            RYi(j,7) = RYi(j,7)*d6;
            RYi(j,8) = RYi(j,8)*d7;
            RYi(j,9) = RYi(j,9)*d8;
            RYi(j,10) = RYi(j,10)*d9;
            RYi(j,11) = RYi(j,11)*d10;
            RYi(j,12) = RYi(j,12)*d11;
            RYi(j,13) = RYi(j,13)*d12;
            IYi(j,1) = IYi(j,1)*d0;
            IYi(j,2) = IYi(j,2)*d1;
            IYi(j,3) = IYi(j,3)*d2;
            IYi(j,4) = IYi(j,4)*d3;
            IYi(j,5) = IYi(j,5)*d4;
            IYi(j,6) = IYi(j,6)*d5;
            IYi(j,7) = IYi(j,7)*d6;
            IYi(j,8) = IYi(j,8)*d7;
            IYi(j,9) = IYi(j,9)*d8;
            IYi(j,10) = IYi(j,10)*d9;
            IYi(j,11) = IYi(j,11)*d10;
            IYi(j,12) = IYi(j,12)*d11;
            IYi(j,13) = IYi(j,13)*d12;
        end
%RYi and IYi are the real and imaginary parts of the spherical harmonics of the
%vectors to all the nearest neighbors of i.  j rows for j nearest neighbors.
%13 columns for -6:m:6
        Rq6mi = (1/NN)*sum(RYi);
        Iq6mi = (1/NN)*sum(IYi);
%Rq6mi and Iq6mi are the sums over all the nearest neighbors, normalized by the
%number of nearest neighbors.  13 rank vectors in m.  q_lm hat in Volkov et al
        magqlm = Rq6mi.*Rq6mi + Iq6mi.*Iq6mi;
        iconst = sqrt(sum(magqlm)); 
%this constant should normalize such that sum_m(q_lm*q_lm) is one
        XRqmIqm(i,1) = xi;
        XRqmIqm(i,2) = yi;
        XRqmIqm(i,3) = zi;
        XRqmIqm(i,4:16) = Rq6mi./iconst;
        XRqmIqm(i,17:29) = Iq6mi./iconst;
        
    end
    radiinew=zeros(numpart,30);
    for i = 1:numpart
        xi = XRqmIqm(i,1);
        yi = XRqmIqm(i,2);
        zi = XRqmIqm(i,3);
        R = sqrt((XRqmIqm(:,1)-xi).*(XRqmIqm(:,1)-xi) + (XRqmIqm(:,2)-yi).*(XRqmIqm(:,2)-yi) + ...
            (XRqmIqm(:,3)-zi).*(XRqmIqm(:,3)-zi));
        radiinew = [XRqmIqm R];
        Neighbornew = radiinew;
        Neighbornew(Neighbornew(:,30)>radius,:)=[];  %includes point i
        [NNN,a]=size(Neighbornew);  %NNN includes point i => number nearest neighbors +1
        DP = zeros(NNN-1,13);
        for j = 1:NNN-1
            for k = 1:13
                DP(j,k)= Neighbornew(1,3+k)*Neighbornew(j+1,3+k)+Neighbornew(1,16+k)*Neighbornew(j+1,16+k);
            end
        end
%DP is the dot product of qlmi and qlmj for every pair of nearest neighbors
%this is now summed over m.
       q6idotq6j = sum(DP,2);
%      connected = q6idotq6j > 0.015; %for fake simple cubic data, 6 nearest neighbors
%      connected = q6idotq6j > 0.035; % for 16 nearest neighbors
%      connected = q6idotq6j > 0.025;  %for 8 nearest neighbors
%      connected = q6idotq6j > 0.03;  %for FCC 12 nearest neighbors
       connected = abs(q6idotq6j) > 0.5;  %for normalized, hopefully
       Q6bonds(i,4) = sum(connected);
 %      Q6allbonds(i,:) = [xi yi zi q6idotq6j'];
end
    fle = ['save ' directory outstem int2str(number) '.dat Q6bonds -ascii'];
    eval(fle);
    number
    %    save file Q6bonds -ascii;
    output = Q6bonds;
 
%output = Q6allbonds;