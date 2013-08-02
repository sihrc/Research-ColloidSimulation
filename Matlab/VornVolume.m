function output=VornVolume(data,scalex,scaley,scalez)

%x y z volume

%faces-edges+vertices=2
%input data has positions in first three columns
   Positions = data(:,1:3);
   Positions(:,1) = Positions(:,1)*scalex;
   Positions(:,2) = Positions(:,2)*scaley;
   Positions(:,3) = Positions(:,3)*scalez;

[nopts,d]=size(Positions);

maxX = max(Positions(:,1));
maxY = max(Positions(:,2));
maxZ = max(Positions(:,3));

[V,C] = voronoin(Positions);
%V is the vertices of the Voronoi diagram, C is a cell array with the
%indices to the positions in V of the vertices that make up the Voronoi cells
%the first point in V is at infinity, if a cell in C contains the index 1,
%it is an unbounded cell

%generates list of Voronoi vertices which lie outside geometrical extent of data set
[lv,wv] = size(V);
badpoints = 1;
for g = 2:lv
    if V(g,1)<0 | V(g,1)>maxX | V(g,2)<0 | V(g,2)>maxY | V(g,3)<0 | V(g,3)>maxZ
        badpoints = [badpoints,g];
    end
end

t = 1;

for i=1:length(C)
    i
    if isempty(intersect(C{i},badpoints))  %get rid of unbounded cells
        matrix = C{i};
        numvert = length(matrix);
        X = V(matrix,1:3);  
        [K,volume] = convhulln(X);
        %K is the indices of the points in X that make up the facets of the cell
        %of course, these are all triangles, with multiple triangles making up a
        %true facet, so now I must write a kluge to put together triangles       
        
        [numtri,d]=size(K);
        
        for z=1:numtri %generate the normal vectors for each triangle
            K(z,4:6)=cross(X(K(z,1),:)-X(K(z,2),:),X(K(z,2),:)-X(K(z,3),:))./...
                norm(cross(X(K(z,1),:)-X(K(z,2),:),X(K(z,2),:)-X(K(z,3),:)));           
        end
        
        faces = zeros(numtri,30);%30 is arbitrary, just a large number of vertices

        for q = 1:numtri
            for r=1:numtri  %find all triangles parallel to the one you're looking at
                K(r,7) = dot(K(q,4:6),K(r,4:6));
            end
            face = K;
            face(:,7) = face(:,7)*10;
            face(:,7) = round(face(:,7));
            face(face(:,7)~=10,:)=[];   
                %this pulls out all triangles with the same normal vector
                %WARNING  this assumes that triangles are drawn in a consistent
                %direction  (ie, signs mean something) THIS MAY NOT BE VALID 
                %there is a check for this written below
            [lef,d]=size(face);
            face = reshape(face(:,1:3),1,3*lef);
            face = unique(face); %now a vector of vertices in actual face, hopefully
            if length(face)>30
                error('face with more than 30 vertices');
            end            
            faces(q,1:length(face)) = face;
        end
        
        faces = unique(faces,'rows'); 
            %faces is a list of all the vertices of the faces, with lots of
            %extra zeros
        %next, a kluge to identify the particle point for the cell in question
        junk = faces;
        [j1,j2]=size(junk);
        junk = reshape(junk,1,j1*j2);
        junk(junk(:)==0)=[];%ooh, this contains the info how many faces share each vertex
        mid1 = X(faces(1,1),:)-(X(faces(1,1),:)-X(faces(1,2),:))./2;
        mid2 = mid1 - (mid1 -X(faces(1,3),:))./2; %a point in the middle of face 1
        notonone = setdiff(junk,faces(1,:)); %list of vertices not on face 1
        midpt = mid2 - (mid2 - X(notonone(1),:))./2;  %should be in cell
        centind = dsearchn(Positions,midpt);
        cent = Positions(centind,:);
        
        res(t,1:3)=cent;
        res(t,4) = volume;
        t = t+1;
    else
        continue
    end           
end 
output = res;
 