function stuff=gofrone(data1); 
%this gofr program actually works.  A little slow, but not too bad.  Eventually
%should adjust this to cross-correlate two different data sets for small/big
%only gives accurate results out to extentX or extentY /2.
	% data1 = dataset number 1, 3+ columns
    % data2 = dataset number 2, 3+ columns
    
    %calculation of g(r) based on description from Chaikin and Lubenski, p37.

    Pos1 = [data1(:,1) data1(:,2) data1(:,3)];
    
    pi = 3.1415926;
    
    maxX1 = max(Pos1(:,1));
    maxY1 = max(Pos1(:,2));
    maxZ1 = max(Pos1(:,3));
    minX1 = min(Pos1(:,1));
    minY1 = min(Pos1(:,2));
    minZ1 = min(Pos1(:,3));
    
    [num1,junk1] = size(Pos1);
    
    vol = (maxX1-minX1)*(maxY1-minY1)*(maxZ1-minZ1); %only approximate for non-square
    dens1 = num1/vol;
    
    r = zeros(num1,1);
dens1
    gofr = zeros(100,4);
    for i = 1:100  %100 steps
        radius = i*0.1; %of 0.1 microns
        inn = radius-0.1;
        count=0;
        vshell = 0;
        for m=1:num1
            if ((Pos1(m,1)-radius < minX1) | (Pos1(m,1)+radius > maxX1))
                continue
            elseif ((Pos1(m,2)-radius < minY1) | (Pos1(m,2)+radius > maxY1))
                continue
            elseif (Pos1(m,3)-radius < minZ1)
                h = Pos1(m,3) - minZ1;
                vshell = vshell + 0.3333333*pi*(2*radius*radius*radius+...
                    3*radius*radius*h) - 0.3333333*pi*(2*inn*inn*inn+3*inn*inn*h);
            elseif (Pos1(m,3)+radius > maxZ1)
                h = maxZ1-Pos1(m,3);
                vshell = vshell + 0.3333333*pi*(2*radius*radius*radius+...
                    3*radius*radius*h) - 0.3333333*pi*(2*inn*inn*inn+3*inn*inn*h);
            else
                vshell = vshell + 4*0.3333*pi*radius*radius*radius - ...
                    4*0.3333*pi*inn*inn*inn;
            end
            for n=1:num1
                r(n) = norm(Pos1(m,:)-Pos1(n,:));
            end
            lessth = [r<radius];
            greaterth = [lessth.*r>inn];
            count = count+sum(greaterth);
        end
        gofr(i,1) = radius;
        gofr(i,2) = count/(vshell*dens1); %maybe off by mult. factor due to dens
        gofr(i,3) = count;
        gofr(i,4) = vshell;
    end
    
 stuff = gofr;
 