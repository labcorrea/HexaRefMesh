%% MeshClass
%
% MeshClass ()
%
% Lucio de Abreu Correa - labcorrea@gmail.com
%
classdef MeshClass < handle
    properties
        % Required Properties
        Coords = [];                % nodes coord  [m]
        NodeStructure = [];         % grid with the nodes in strucutred topology
        Nel =[];                    % vector with the number of elements in [x y z]
        Connectivity = [];          % connectivity of the elements
        ElStructure = [];           % grid with the elements in strucutred topology
        Father = [];                % the name of the father element in the ref procedure
    end
    methods
        %% CONSTRUCTOR
        function obj = MeshClass( )
            %% MeshClass
            % MeshClass contructor
            %
            % Syntax:
            %   newobj = MeshClass (  );
            %
            % Inputs:
            %
            % See also
        end
        function newobj = copyobj(obj)
            %% copyobj
            % Return a new class instance containing the same properties
            % values from the input
            %
            % Syntax
            %   newobj = copyobj(obj)
            %
            % Inputs:
            %   obj : Object to be copied
            %
            % Outputs:
            %   newobj : A new object from the same class with a copy of
            %   all its properties
            %
            if isscalar(obj)
                newobj = eval(mfilename);
                props = properties(obj);
                for i_props = 1:numel(props)
                    if isobject(obj.(props{i_props}))
                        if ~isempty(obj.(props{i_props}))
                            newobj.(props{i_props}) = obj.(props{i_props}).copyobj;
                        else
                            %empty
                        end
                    else
                        newobj.(props{i_props}) = obj.(props{i_props});
                    end
                end
            else
                for i=1:numel(obj)
                    newobj(i) = obj(i).copyobj; %#ok<AGROW>
                end
                newobj = reshape(newobj,size(obj));
            end
        end
        %% BUILD MESH
        function obj = Mesher(obj,Geometry)
            %% MeshClass.Mesher
            % Build strucutred mesh for plates and pipes
            %
            % Syntax:
            %   newobj = obj.mesher;
            %
            % Inputs:
            %   Geometry.Dimension = [x y z]
            % Outputs:
            %
            % See also,
            %
            if prod(obj.Nel) > 2^32
                error('Please reduce the number of elements')
            end

            %nodes
            nnodes = obj.Nel + 1;
            obj.NodeStructure = zeros(prod(double(nnodes)),4);
            index = 1;
            for k = 1 : nnodes(3)
                for j = 1 : nnodes(2)
                    for i = 1 : nnodes(1)
                        name = i + (j-1)*nnodes(1) + (k-1)*nnodes(2)*nnodes(1);
                        obj.NodeStructure(index,:) = [name i j k];
                        index = index + 1 ;
                    end
                end
            end
            clear index name
            %elements
            index = 1;
            obj.Connectivity = zeros(prod(double(obj.Nel)),9);
            obj.Father =  zeros(prod(double(obj.Nel)),1);
            for k = 1 : obj.Nel(3)
                for j = 1 : obj.Nel(2)
                    for i = 1 :obj.Nel(1)
                        ElId = i + (j-1)*obj.Nel(1) + (k-1)*obj.Nel(2)*obj.Nel(1);
                        obj.Connectivity(index,:) = [ElId ...
                            (i+0) + (j+0-1)*nnodes(1) + (k+0-1)*nnodes(2)*nnodes(1) ...
                            (i+1) + (j+0-1)*nnodes(1) + (k+0-1)*nnodes(2)*nnodes(1) ...
                            (i+1) + (j+1-1)*nnodes(1) + (k+0-1)*nnodes(2)*nnodes(1) ...
                            (i+0) + (j+1-1)*nnodes(1) + (k+0-1)*nnodes(2)*nnodes(1) ...
                            (i+0) + (j+0-1)*nnodes(1) + (k+1-1)*nnodes(2)*nnodes(1) ...
                            (i+1) + (j+0-1)*nnodes(1) + (k+1-1)*nnodes(2)*nnodes(1) ...
                            (i+1) + (j+1-1)*nnodes(1) + (k+1-1)*nnodes(2)*nnodes(1) ...
                            (i+0) + (j+1-1)*nnodes(1) + (k+1-1)*nnodes(2)*nnodes(1)];
                        obj.ElStructure(index,:) = [ ElId, i, j, k];
                        obj.Father(index) = ElId;
                        index = index + 1 ;
                    end
                end
                if(mod(100*double(k)/double(obj.Nel(3)),10)==0)
                    fprintf('Building Connectivity, %f complete\n',100*double(k)/double(obj.Nel(3)))
                end
            end
            %coords
            x = Geometry.Dimension(1)*double(obj.NodeStructure(:,2)-1)/double(max(obj.NodeStructure(:,2))-1);
            y = Geometry.Dimension(2)*double(obj.NodeStructure(:,3)-1)/double(max(obj.NodeStructure(:,3))-1);
            z = Geometry.Dimension(3)*double(obj.NodeStructure(:,4)-1)/double(max(obj.NodeStructure(:,4))-1);

            obj.Coords = [x y z];
        end
        %% NODE functions
        function name = AddNode(obj,node, coord)
            %% MeshClass.AddNode
            % create a new nodes
            %
            % Syntax:
            %   obj = obj.AddNode(node, coord);
            %
            % Inputs:
            %   node : vector 1 by 3 with ref coords
            %   coord: vector 1 by 3 with coords
            %
            % Outputs:
            %
            % See also,
            %
            aux = logical((obj.NodeStructure(:,2) == node(1)) .* ...
                (obj.NodeStructure(:,3) == node(2)) .* ...
                (obj.NodeStructure(:,4) == node(3)));

            switch sum(aux)
                case 0
                    %create the node
                    obj.NodeStructure(end+1,:) = [(size(obj.NodeStructure,1)+1) node];
                    obj.Coords(end+1,:) = coord;
                    name = size(obj.Coords,1);
                case 1
                    %node already exist
                    name = obj.NodeStructure(aux,1);
                otherwise
                    error('The sum gave a not possible number :)');
            end
        end
        %% ELEMENT functions
        function ElId = GetElementNeighbors(obj,iel,mode)
            %% MeshClass.GetElementNeighbors
            % Get Neighbors Elements
            %
            % Syntax:
            %   obj = obj.GetElementNeighbors(iel, mode);
            %
            % Inputs:
            %   iel : element id
            %   mode: in the plane or all
            %
            % Outputs:
            %
            % See also,
            %
            elref = obj.ElStructure(iel,2:end);
            x(1) = elref(1)-3;
            x(2) = elref(1)+3;
            y(1) = elref(2)-3;
            y(2) = elref(2)+3;
            z(1) = elref(3)-3;
            z(2) = elref(3)+3;
            xx = (obj.ElStructure(:,2) >= x(1)) .* (obj.ElStructure(:,2) <= x(2));
            yy = (obj.ElStructure(:,3) >= y(1)) .* (obj.ElStructure(:,3) <= y(2));
            switch mode
                case 'all'
                    zz = (obj.ElStructure(:,4) >= z(1)) .* (obj.ElStructure(:,4) <= z(2));
                case 'plane'
                    zz = obj.ElStructure(:,4) == elref(3);
            end
            ind = logical(xx.*yy.*zz);
            ElId = obj.ElStructure(ind,1);
        end
        %% MESH Refinement
        function [NList,ElList] = BB(obj,BoundBox)
            %% MeshClass.BB
            % Found the bounding box and return all the elements inside of
            % the BB
            %
            % Syntax:
            %   obj = obj.BB( BoundBox );
            %
            % Inputs:
            %   BoundBox : BB with [xmin xmax; ymin ymax;zmin zmax]
            %
            % Outputs:
            %   ElList : Logical list of the elements touched by the BB
            %   NList : Logical list of the nodes inside the BB
            %
            % See also,
            NList = logical((obj.Coords(:,1)>=BoundBox(1,1)).*(obj.Coords(:,1)<=BoundBox(1,2)) ...
                .*(obj.Coords(:,2)>=BoundBox(2,1)).*(obj.Coords(:,2)<=BoundBox(2,2)) ...
                .*(obj.Coords(:,3)>=BoundBox(3,1)).*(obj.Coords(:,3)<=BoundBox(3,2)));
            if nargout == 2
                nodeName= obj.NodeStructure(NList,1);
                aux = obj.Connectivity(:,2:end);
                aux = reshape(aux',1, 8*size(obj.Connectivity,1))';
                El = zeros(numel(aux),1);
                for inode = 1:numel(nodeName)
                    El = (aux == nodeName(inode)) + El;
                end
                ElList = El>0;
                clear El
                ElList = reshape(ElList,[8 size(obj.Connectivity,1)]);
                ElList = sum(ElList,1)>0;
                ElList = ElList.';
            end
        end
        function RefMesh(obj,BoundBox)
            %% MeshClass.RefMesh
            % Apply a 27 refinement in the BB region and make spread the
            % changes to all keep the conform hexaedrical mesh.
            %
            % Syntax:
            %   obj = obj.RefMesh( BoundBox );
            %
            % Inputs:
            %   BoundBox : BB with [xmin xmax; ymin ymax;zmin zmax]
            %
            % Outputs:
            %
            % See also,

            %% update the RefNodes and RefElements
            obj.NodeStructure(:,2:end) = 3*(obj.NodeStructure(:,2:end)-1)+1;
            obj.ElStructure(:,2:end) = 3*(obj.ElStructure(:,2:end)-1)+1;
            %% found the bouding box
            [~,ElList] = obj.BB(BoundBox);
            %% create the new 27-tree
            refEl = [];
            for iel = 1: numel(ElList)
                if(ElList(iel))
                    obj.HexaTo27Tree(iel);
                    refEl = [refEl iel];
                end
            end
            if isscalar(refEl)
                error('It should have more than one element');
            end
            %% Make the mesh conform
            di = [];
            for iel = 1 : numel(refEl)
                ids = obj.GetElementNeighbors(refEl(iel),'plane');
                di = [di setdiff(ids,refEl)'];
            end

            %split the face and edge ref
            diUnique = unique(di);
            faceRef = [];
            edgeRef = [];
            for iel = 1 : numel(diUnique)
                switch sum(di == diUnique(iel))
                    case 1
                        edgeRef = [edgeRef diUnique(iel)];
                    otherwise
                        faceRef = [faceRef diUnique(iel)];
                end
            end

            for iel = 1 : numel(faceRef)
                ids = obj.GetElementNeighbors(faceRef(iel),'plane');
                a = obj.ElStructure(faceRef(iel),:);
                b = obj.ElStructure(setdiff(ids,faceRef(iel)),:);

                c = 1;
                for jel = 1:size(b,1)
                    if sum(refEl==b(jel,1))~=0
                        aux(c,:) = b(jel,2:end) - a(2:end);
                        c = c+1;
                    end
                end

                aux = sum(aux);
                [a ,b] = max(abs(aux));
                a = aux(b);
                if sign(a)==-1
                    if b == 2
                        iface = 3;
                    elseif b == 1
                        iface = 6;
                    else
                        error('Not possible')
                    end
                end
                if sign(a)==1
                    if b == 2
                        iface = 5;
                    elseif b == 1
                        iface = 4;
                    else
                        error('Not possible')
                    end
                end
                obj.HexaTo27TreeTemplateFace(faceRef(iel),iface);
            end

            for iel = 1 : numel(edgeRef)
                ids = obj.GetElementNeighbors(edgeRef(iel),'plane');
                a = obj.ElStructure(edgeRef(iel),:);
                b = obj.ElStructure(setdiff(ids,edgeRef(iel)),:);

                c = 1;
                for jel = 1:size(b,1)
                    if sum(refEl==b(jel,1))~=0
                        aux(c,:) = b(jel,2:end) - a(2:end);
                        c = c+1;
                    end
                end

                if size(aux,1)~=1
                    error('Corner Element with more than one... something is wrong =0')
                end
                if sign(aux(1))==1 && sign(aux(2))==1
                    iedge = 7;
                end
                if sign(aux(1))==-1 && sign(aux(2))==-1
                    iedge = 5;
                end
                if sign(aux(1))==1 && sign(aux(2))==-1
                    iedge = 6;
                end
                if sign(aux(1))==-1 && sign(aux(2))==1
                    iedge = 8;
                end

                obj.HexaTo27TreeTemplateEdge(edgeRef(iel),iedge);
            end

        end
        function HexaTo27Tree(obj,iel)
            %% MeshClass.HexaTo27Tree
            % Apply a 27 refinement one element
            %
            % Syntax:
            %   obj = obj.HexaTo27Tree( iel );
            %
            % Inputs:
            %   iel : element name
            %
            % Outputs:
            %
            % See also,

            %build the template to the elements
            step = 3;

            %reference element edge 0
            local_ref = zeros(27,8,3);
            count = 1;
            for iz = [0 3 6]
                for iy = [0 3 6]
                    for ix = [0 3 6]

                        x = ix;
                        y = iy;
                        z = iz;

                        local_ref(count,1,:) = [x      y      z];
                        local_ref(count,2,:) = [x+step y      z];
                        local_ref(count,3,:) = [x+step y+step z];
                        local_ref(count,4,:) = [x      y+step z];

                        local_ref(count,5,:) = [x      y      z+step];
                        local_ref(count,6,:) = [x+step y      z+step];
                        local_ref(count,7,:) = [x+step y+step z+step];
                        local_ref(count,8,:) = [x      y+step z+step];

                        count = count + 1;
                    end
                end
            end
            local_ref = 2*local_ref/9;
            local_ref = local_ref -1;

            ord = 1:8;
            %find the element and nodes
            ElMask = obj.Connectivity(:,1) == iel;
            nodes = obj.Connectivity(ElMask,2:9);
            aux = zeros(8,4);
            for ino = 1 : 8
                aux(ino,:) = obj.NodeStructure(nodes(ino) ,:);
            end
            rotatedNodes = nodes(ord);
            %get the coords of the nodes
            coords = obj.Coords(rotatedNodes,:);

            index = size(obj.Connectivity,1);
            nodeName = zeros(1,8);
            base = [min(aux(:,2)) min(aux(:,3)) min(aux(:,4))];

            for ie = 1 : 27
                %build a ref coords for a [-1 1] cubic element
                ref = squeeze(local_ref(ie,:,:));
                refno = floor((ref+1)*3/2) + base;

                % make the interpolation
                MapCoord = obj.LinearMap(ref,coords);
                for inode = 1 : 8
                    nodeName(inode) = obj.AddNode(refno(inode,:),MapCoord(inode,:));
                end
                if ie == 1
                    nameEl = iel;
                else
                    index = index + 1 ;
                    nameEl = index;
                end
                obj.Connectivity(nameEl,:) = [nameEl nodeName];
                obj.Father(nameEl) = iel;
            end
        end
        function HexaTo27TreeTemplateFace(obj,iel,iFace)
            %% MeshClass.HexaTo27TreeTemplateEdge
            % Make the element iel a face template to avoid hanging nodes
            %
            % Syntax:
            %   obj = obj.HexaTo27TreeTemplateface( iel, face );
            %
            % Inputs:
            %   iel  : element name
            %   iFace: face to be refined
            %
            % Outputs:
            %
            % See also,

            %build the template to the elements
            step = 2/3;

            %reference element edge 0
            local_ref = zeros(13,8,3);
            %element 0
            inel = 1;
            local_ref(inel,1,:) = [-1 -1 -1];
            local_ref(inel,2,:) = [ 1 -1 -1];
            local_ref(inel,3,:) = [ 1  1 -1];
            local_ref(inel,4,:) = [-1  1 -1];

            local_ref(inel,5,:) = [-1 -1+step  -1+step];
            local_ref(inel,6,:) = [ 1 -1+step  -1+step];
            local_ref(inel,7,:) = [ 1  -1+2*step  -1+step];
            local_ref(inel,8,:) = [-1  -1+2*step  -1+step];
            %element 1
            inel = 2;
            local_ref(inel,1,:) = [-1 -1+step -1+step];
            local_ref(inel,2,:) = [ 1 -1+step -1+step];
            local_ref(inel,3,:) = [ 1  -1+2*step -1+step];
            local_ref(inel,4,:) = [-1  -1+2*step -1+step];

            local_ref(inel,5,:) = [-1+step -1+step  -1+2*step];
            local_ref(inel,6,:) = [ -1+2*step -1+step  -1+2*step];
            local_ref(inel,7,:) = [ -1+2*step  -1+2*step  -1+2*step];
            local_ref(inel,8,:) = [-1+step  -1+2*step  -1+2*step];
            %element 2
            inel = 3;
            local_ref(inel,1,:) = [-1+step -1+step -1+2*step];
            local_ref(inel,2,:) = [-1+2*step -1+step -1+2*step];
            local_ref(inel,3,:) = [-1+2*step -1+2*step -1+2*step];
            local_ref(inel,4,:) = [-1+step -1+2*step -1+2*step];

            local_ref(inel,5,:) = [-1+step -1+step 1];
            local_ref(inel,6,:) = [-1+2*step -1+step 1];
            local_ref(inel,7,:) = [-1+2*step -1+2*step 1];
            local_ref(inel,8,:) = [-1+step -1+2*step 1];
            %element 3
            inel = 4;
            local_ref(inel,1,:) = [-1  -1  -1];
            local_ref(inel,2,:) = [1 -1 -1];
            local_ref(inel,3,:) = [1  -1+step  -1+step];
            local_ref(inel,4,:) = [-1  -1+step -1+step];

            local_ref(inel,5,:) = [-1+step -1  -1+step];
            local_ref(inel,6,:) = [-1+2*step  -1 -1+step];
            local_ref(inel,7,:) = [-1+2*step -1+1*step  -1+2*step];
            local_ref(inel,8,:) = [-1+step  -1+1*step -1+2*step];
            %element 4
            inel = 5;
            local_ref(inel,1,:) = [-1  -1+2*step  -1+step];
            local_ref(inel,2,:) = [1  -1+2*step  -1+step];
            local_ref(inel,3,:) = [1   1 -1];
            local_ref(inel,4,:) = [-1  1 -1];

            local_ref(inel,5,:) = [-1+step  -1+2*step  -1+2*step];
            local_ref(inel,6,:) = [-1+2*step  -1+2*step  -1+2*step];
            local_ref(inel,7,:) = [-1+2*step 1  -1+step];
            local_ref(inel,8,:) = [-1+step 1  -1+step];
            %element 5
            inel = 6;
            local_ref(inel,1,:) = [-1 -1 -1];
            local_ref(inel,2,:) = [-1+step -1 -1+step];
            local_ref(inel,3,:) = [-1+step -1+step -1+2*step];
            local_ref(inel,4,:) = [-1 -1+step -1+step];

            local_ref(inel,5,:) = [-1 -1 1];
            local_ref(inel,6,:) = [-1+step -1 1];
            local_ref(inel,7,:) = [-1+step -1+step 1];
            local_ref(inel,8,:) = [-1 -1+step 1];
            %element 6
            inel = 7;
            local_ref(inel,1,:) = [-1+step -1 -1+1*step];
            local_ref(inel,2,:) = [-1+2*step -1 -1+step];
            local_ref(inel,3,:) = [-1+2*step -1+step -1+2*step];
            local_ref(inel,4,:) = [-1+step -1+step -1+2*step];

            local_ref(inel,5,:) = [-1+step -1 1];
            local_ref(inel,6,:) = [-1+2*step -1 1];
            local_ref(inel,7,:) = [-1+2*step -1+1*step 1];
            local_ref(inel,8,:) = [-1+step -1+1*step 1];
            %element 7
            inel = 8;
            local_ref(inel,1,:) = [-1+2*step -1 -1+step];
            local_ref(inel,2,:) = [1 -1 -1];
            local_ref(inel,3,:) = [1 -1+step -1+1*step];
            local_ref(inel,4,:) = [-1+2*step -1+1*step -1+2*step];

            local_ref(inel,5,:) = [-1+2*step -1 1];
            local_ref(inel,6,:) = [1 -1 1];
            local_ref(inel,7,:) = [1 -1+1*step 1];
            local_ref(inel,8,:) = [-1+2*step -1+1*step 1];
            %element 8
            inel = 9;
            local_ref(inel,1,:) = [-1+2*step -1+1*step -1+2*step];
            local_ref(inel,2,:) = [1 -1+1*step -1+1*step];
            local_ref(inel,3,:) = [1 -1+2*step -1+1*step];
            local_ref(inel,4,:) = [-1+2*step -1+2*step -1+2*step];

            local_ref(inel,5,:) = [-1+2*step -1+1*step 1];
            local_ref(inel,6,:) = [1 -1+1*step 1];
            local_ref(inel,7,:) = [1 -1+2*step 1];
            local_ref(inel,8,:) = [-1+2*step -1+2*step 1];
            %element 9
            inel = 10;
            local_ref(inel,1,:) = [-1+2*step -1+2*step -1+2*step];
            local_ref(inel,2,:) = [1 -1+2*step -1+1*step];
            local_ref(inel,3,:) = [1  1 -1];
            local_ref(inel,4,:) = [-1+2*step  1 -1+1*step];

            local_ref(inel,5,:) = [-1+2*step -1+2*step 1];
            local_ref(inel,6,:) = [1 -1+2*step 1];
            local_ref(inel,7,:) = [1 1 1];
            local_ref(inel,8,:) = [-1+2*step 1 1];
            %element 10
            inel = 11;
            local_ref(inel,1,:) = [-1+1*step -1+2*step -1+2*step];
            local_ref(inel,2,:) = [-1+2*step -1+2*step -1+2*step];
            local_ref(inel,3,:) = [-1+2*step 1 -1+step];
            local_ref(inel,4,:) = [-1+1*step  1 -1+1*step];

            local_ref(inel,5,:) = [-1+1*step -1+2*step 1];
            local_ref(inel,6,:) = [-1+2*step -1+2*step 1];
            local_ref(inel,7,:) = [-1+2*step 1 1];
            local_ref(inel,8,:) = [-1+1*step 1 1];
            %element 11
            inel = 12;
            local_ref(inel,1,:) = [-1 -1+2*step -1+1*step];
            local_ref(inel,2,:) = [-1+1*step -1+2*step -1+2*step];
            local_ref(inel,3,:) = [-1+1*step  1 -1+step];
            local_ref(inel,4,:) = [-1  1 -1];

            local_ref(inel,5,:) = [-1 -1+2*step 1];
            local_ref(inel,6,:) = [-1+1*step -1+2*step 1];
            local_ref(inel,7,:) = [-1+1*step 1 1];
            local_ref(inel,8,:) = [-1 1 1];
            %element 12
            inel = 13;
            local_ref(inel,1,:) = [-1 -1+1*step -1+1*step];
            local_ref(inel,2,:) = [-1+1*step -1+1*step -1+2*step];
            local_ref(inel,3,:) = [-1+1*step -1+2*step -1+2*step];
            local_ref(inel,4,:) = [-1 -1+2*step -1+step];

            local_ref(inel,5,:) = [-1 -1+1*step 1];
            local_ref(inel,6,:) = [-1+1*step -1+1*step 1];
            local_ref(inel,7,:) = [-1+1*step -1+2*step 1];
            local_ref(inel,8,:) = [-1 -1+2*step 1];

            switch iFace
                case 1
                    %edge 0 1 2 3
                    rot = [0 0 0];
                    sym = [0 0 1];
                    rotOrd = [1 2 3];
                case 2
                    %edge 8 9 10 11
                    rot = [0 0 0];
                    sym = [0 0 0];
                    rotOrd = [1 2 3];
                case 3
                    %edge 0 4 5 8
                    rot = [-1 0 0];
                    sym = [0 0 0];
                    rotOrd = [1 3 2];
                case 4
                    %edge 1 5 6 9
                    rot = [0 -1 0];
                    sym = [1 0 0];
                    rotOrd = [3 2 1];
                case 5
                    %edge 2 6 7 10
                    rot = [1 0 0];
                    sym = [0 1 0];
                    rotOrd = [1 3 2];
                case 6
                    %edge 3 4 7 11
                    rot = [0 1 0];
                    sym = [0 0 0];
                    rotOrd = [3 2 1];
            end
            ord = MeshClass.RotateHex(rot,sym);

            %find the element and nodes
            ElMask = obj.Connectivity(:,1) == iel;
            nodes = obj.Connectivity(ElMask,2:9);
            aux = zeros(8,4);
            for ino = 1 : 8
                aux(ino,:) = obj.NodeStructure(nodes(ino) == obj.NodeStructure(:,1),:);
            end
            rotatedNodes = nodes(ord);
            %get the coords of the nodes
            coords = obj.Coords(rotatedNodes,:);
            index = size(obj.Connectivity,1);
            for ie = 1:13

                %build a ref coords for a [-1 1] cubic element
                ref = (squeeze(local_ref(ie,:,:)));
                % make the interpolation
                MapCoord = obj.LinearMap(ref,coords);

                base = [min(aux(ord,2)) min(aux(ord,3)) min(aux(ord,4))];
                refno = floor((ref(:,:)+1)*3/2) + base(rotOrd);
                nodeName = zeros(1,8);
                for inode = 1 : 8
                    if((ref(inode,1)==1 || ref(inode,1)==-1) && (ref(inode,2)==1 || ref(inode,2)==-1) &&(ref(inode,3)==1 || ref(inode,3)==-1))
                        if (ref(inode,1)==-1) && (ref(inode,2)==-1) && (ref(inode,3)==-1)
                            nodeName(inode)= rotatedNodes(1);
                        elseif (ref(inode,1)==1) && (ref(inode,2)==-1) && (ref(inode,3)==-1)
                            nodeName(inode)= rotatedNodes(2);
                        elseif (ref(inode,1)==1) && (ref(inode,2)==1) && (ref(inode,3)==-1)
                            nodeName(inode)= rotatedNodes(3);
                        elseif (ref(inode,1)==-1) && (ref(inode,2)==1) && (ref(inode,3)==-1)
                            nodeName(inode)= rotatedNodes(4);
                        elseif (ref(inode,1)==-1) && (ref(inode,2)==-1) && (ref(inode,3)==1)
                            nodeName(inode)= rotatedNodes(5);
                        elseif (ref(inode,1)==1) && (ref(inode,2)==-1) && (ref(inode,3)==1)
                            nodeName(inode)= rotatedNodes(6);
                        elseif (ref(inode,1)==1) && (ref(inode,2)==1) && (ref(inode,3)==1)
                            nodeName(inode)= rotatedNodes(7);
                        elseif (ref(inode,1)==-1) && (ref(inode,2)==1) && (ref(inode,3)==1)
                            nodeName(inode)= rotatedNodes(8);
                        end
                    else
                        nodeName(inode) = obj.AddNode(refno(inode,rotOrd),MapCoord(inode,:));
                    end
                end
                if ie == 1
                    nameEl = iel;
                else
                    index = index + 1 ;
                    nameEl = index;
                end
                obj.Connectivity(nameEl,:) = [nameEl nodeName(ord)];
                obj.Father(nameEl) = iel;
            end
        end
        function HexaTo27TreeTemplateEdge(obj,iel,iEdge)
            %% MeshClass.HexaTo27TreeTemplateEdge
            % Make the element iel a edge template to avoid hanging nodes
            %
            % Syntax:
            %   obj = obj.HexaTo27TreeTemplateEdge( iel, edge );
            %
            % Inputs:
            %   iel  : element name
            %   iEdge: edge to be refined
            %
            % Outputs:
            %
            % See also,
            step = 2/3;
            %reference element edge 0
            local_ref = zeros(5,8,3);
            %element 0
            inel = 1;
            local_ref(inel,1,:) = [-1 -1 -1];
            local_ref(inel,2,:) = [-1+step -1 -1];
            local_ref(inel,3,:) = [-1+step -1+2*step -1];
            local_ref(inel,4,:) = [-1 1 -1];

            local_ref(inel,5,:) = [-1 -1 1];
            local_ref(inel,6,:) = [-1+step -1 -1+2*step];
            local_ref(inel,7,:) = [-1+step -1+2*step -1+2*step];
            local_ref(inel,8,:) = [-1 1 1];
            %element 1
            inel = 2;
            local_ref(inel,1,:) = [-1+step -1 -1];
            local_ref(inel,2,:) = [-1+2*step -1 -1];
            local_ref(inel,3,:) = [-1+2*step -1+2*step -1];
            local_ref(inel,4,:) = [-1+step -1+2*step -1];

            local_ref(inel,5,:) = [-1+step -1 -1+2*step];
            local_ref(inel,6,:) = [-1+2*step -1 -1+2*step];
            local_ref(inel,7,:) = [-1+2*step -1+2*step -1+2*step];
            local_ref(inel,8,:) = [-1+step -1+2*step -1+2*step];
            %element 2
            inel = 3;
            local_ref(inel,1,:) = [-1+2*step -1 -1];
            local_ref(inel,2,:) = [1 -1 -1];
            local_ref(inel,3,:) = [1  1 -1];
            local_ref(inel,4,:) = [-1+2*step -1+2*step -1];

            local_ref(inel,5,:) = [-1+2*step -1 -1+2*step];
            local_ref(inel,6,:) = [1 -1  1];
            local_ref(inel,7,:) = [1 1 1];
            local_ref(inel,8,:) = [-1+2*step -1+2*step -1+2*step];
            %element 3
            inel = 4;
            local_ref(inel,1,:) = [-1+1*step -1+2*step -1];
            local_ref(inel,2,:) = [-1+2*step -1+2*step -1];
            local_ref(inel,3,:) = [1  1 -1];
            local_ref(inel,4,:) = [-1 1 -1];

            local_ref(inel,5,:) = [-1+1*step -1+2*step -1+2*step];
            local_ref(inel,6,:) = [-1+2*step -1+2*step -1+2*step];
            local_ref(inel,7,:) = [1 1 1];
            local_ref(inel,8,:) = [-1 1 1];
            %element 4
            inel = 5;
            local_ref(inel,1,:) = [-1+step -1 -1+2*step];
            local_ref(inel,2,:) = [-1+2*step -1 -1+2*step];
            local_ref(inel,3,:) = [-1+2*step -1+2*step -1+2*step];
            local_ref(inel,4,:) = [-1+step -1+2*step -1+2*step];

            local_ref(inel,5,:) = [-1 -1 1];
            local_ref(inel,6,:) = [1 -1 1];
            local_ref(inel,7,:) = [1 1 1];
            local_ref(inel,8,:) = [-1 1 1];

            switch iEdge
                %                 case 1
                %                     rot = [0 0 0];
                %                     sym = [0 0 0];
                %                     rotOrd = [1 2 3];
                %                 case 2
                %                     rot = [0 0 -1];
                %                     sym = [0 0 0];
                %                     rotOrd = [1 2 3];
                %                 case 3
                %                     rot = [0 0 0];
                %                     sym = [0 1 0];
                %                     rotOrd = [1 2 3];
                %                 case 4
                %                     rot = [0 0 1];
                %                     sym = [0 0 0];
                %                     rotOrd = [2 1 3];
                case 5
                    rot = [0 -1 0];
                    sym = [0 0 0];
                    rotOrd = [1 2 3];
                case 6
                    rot = [0 1 0];
                    sym = [0 0 0];
                    rotOrd = [2 1 3];
                case 7
                    rot = [0 1 0];
                    sym = [0 1 0];
                    rotOrd = [2 1 3];
                case 8
                    rot = [0 -1 0];
                    sym = [0 1 0];
                    rotOrd = [1 2 3];
                    %                 case 9
                    %                     rot = [0 0 0];
                    %                     sym = [0 0 1];
                    %                     rotOrd = [1 2 3];
                    %                 case 10
                    %                     rot = [0 0 -1];
                    %                     sym = [0 0 1];
                    %                     rotOrd = [1 2 3];
                    %                 case 11
                    %                     rot = [-1 0 0];
                    %                     sym = [0 1 0];
                    %                     rotOrd = [1 2 3];
                    %                 case 12
                    %                     rot = [0 0 1];
                    %                     sym = [0 0 1];
                    %                     rotOrd = [1 2 3];
            end
            ord = MeshClass.RotateHex(rot,sym);

            %find the element and nodes
            ElMask = obj.Connectivity(:,1) == iel;
            nodes = obj.Connectivity(ElMask,2:9);
            aux = zeros(8,4);
            for ino = 1 : 8
                aux(ino,:) = obj.NodeStructure(nodes(ino) == obj.NodeStructure(:,1),:);
            end
            rotatedNodes = nodes(ord);
            %get the coords of the nodes
            coords = obj.Coords(rotatedNodes,:);
            index = size(obj.Connectivity,1);
            for ie = 1:5

                %build a ref coords for a [-1 1] cubic element
                ref = (squeeze(local_ref(ie,:,:)));
                % make the interpolation
                MapCoord = obj.LinearMap(ref,coords);

                base = [min(aux(ord,2)) min(aux(ord,3)) min(aux(ord,4))];
                refno = floor((ref(:,:)+1)*3/2) + base(rotOrd);
                nodeName = zeros(1,8);
                for inode = 1 : 8
                    if((ref(inode,1)==1 || ref(inode,1)==-1) && (ref(inode,2)==1 || ref(inode,2)==-1) &&(ref(inode,3)==1 || ref(inode,3)==-1))
                        if (ref(inode,1)==-1) && (ref(inode,2)==-1) && (ref(inode,3)==-1)
                            nodeName(inode)= rotatedNodes(1);
                        elseif (ref(inode,1)==1) && (ref(inode,2)==-1) && (ref(inode,3)==-1)
                            nodeName(inode)= rotatedNodes(2);
                        elseif (ref(inode,1)==1) && (ref(inode,2)==1) && (ref(inode,3)==-1)
                            nodeName(inode)= rotatedNodes(3);
                        elseif (ref(inode,1)==-1) && (ref(inode,2)==1) && (ref(inode,3)==-1)
                            nodeName(inode)= rotatedNodes(4);
                        elseif (ref(inode,1)==-1) && (ref(inode,2)==-1) && (ref(inode,3)==1)
                            nodeName(inode)= rotatedNodes(5);
                        elseif (ref(inode,1)==1) && (ref(inode,2)==-1) && (ref(inode,3)==1)
                            nodeName(inode)= rotatedNodes(6);
                        elseif (ref(inode,1)==1) && (ref(inode,2)==1) && (ref(inode,3)==1)
                            nodeName(inode)= rotatedNodes(7);
                        elseif (ref(inode,1)==-1) && (ref(inode,2)==1) && (ref(inode,3)==1)
                            nodeName(inode)= rotatedNodes(8);
                        end
                    else
                        nodeName(inode) = obj.AddNode(refno(inode,rotOrd),MapCoord(inode,:));
                    end
                end
                if ie == 1
                    nameEl = iel;
                else
                    index = index + 1 ;
                    nameEl = index;
                end
                obj.Connectivity(nameEl,:) = [nameEl nodeName(ord)];
                obj.Father(nameEl) = iel;
            end
        end
        function MapCoord = LinearMap(obj,ref,coords)
            %% MeshClass.LinearMap
            % Provide a linear map/interpolation in a reference element
            %
            % Syntax:
            %   obj = LinearMap(obj,ref,coords)
            %
            % Inputs:
            %   ref : vector N by 3 ref in a [-1 1] cube
            %   coords: vector 8x3 coords of the exterior hexaedra
            %
            % Outputs:
            %   MapCoord : vector of N by 3 coords
            %
            % See also,

            N(:,1) = (1-ref(:,1)).*(1-ref(:,2)).*(1-ref(:,3))*0.125;
            N(:,2) = (1+ref(:,1)).*(1-ref(:,2)).*(1-ref(:,3))*0.125;
            N(:,3) = (1+ref(:,1)).*(1+ref(:,2)).*(1-ref(:,3))*0.125;
            N(:,4) = (1-ref(:,1)).*(1+ref(:,2)).*(1-ref(:,3))*0.125;

            N(:,5) = (1-ref(:,1)).*(1-ref(:,2)).*(1+ref(:,3))*0.125;
            N(:,6) = (1+ref(:,1)).*(1-ref(:,2)).*(1+ref(:,3))*0.125;
            N(:,7) = (1+ref(:,1)).*(1+ref(:,2)).*(1+ref(:,3))*0.125;
            N(:,8) = (1-ref(:,1)).*(1+ref(:,2)).*(1+ref(:,3))*0.125;

            MapCoord = N*coords;
        end
        %% plot function
        function MeshScatter(obj,h)

            if exist('h','var')
                figure(h)
            else
                h = figure;
            end
            switch size(obj.Coords,2)
                case 1
                    plot(ones(numel(obj.Coords(:,1)),1),obj.Coords(:,1),'x-');
                case 2
                    scatter(obj.Coords(:,1),obj.Coords(:,2))
                case 3
                    scatter3(obj.Coords(:,1),obj.Coords(:,2),obj.Coords(:,3))
                otherwise
                    error('Dim not supported')
            end

        end
        function PlotMesh(obj,h)
            if exist('h','var')
                figure(h)
            else
                h = figure;
            end
            for iel = 1 : size(obj.Connectivity,1)
                switch size(obj.Coords,2)
                    case 3
                        nodes = obj.Connectivity(iel,2:9);
                        coords = obj.Coords(nodes,:);

                        line([coords(1,1) coords(2,1)],[coords(1,2) coords(2,2)],[coords(1,3) coords(2,3)],'Color','black')
                        line([coords(2,1) coords(3,1)],[coords(2,2) coords(3,2)],[coords(2,3) coords(3,3)],'Color','black')
                        line([coords(3,1) coords(4,1)],[coords(3,2) coords(4,2)],[coords(3,3) coords(4,3)],'Color','black')
                        line([coords(4,1) coords(1,1)],[coords(4,2) coords(1,2)],[coords(4,3) coords(1,3)],'Color','black')

                        line([coords(5,1) coords(6,1)],[coords(5,2) coords(6,2)],[coords(5,3) coords(6,3)],'Color','black')
                        line([coords(6,1) coords(7,1)],[coords(6,2) coords(7,2)],[coords(6,3) coords(7,3)],'Color','black')
                        line([coords(7,1) coords(8,1)],[coords(7,2) coords(8,2)],[coords(7,3) coords(8,3)],'Color','black')
                        line([coords(8,1) coords(5,1)],[coords(8,2) coords(5,2)],[coords(8,3) coords(5,3)],'Color','black')

                        line([coords(1,1) coords(5,1)],[coords(1,2) coords(5,2)],[coords(1,3) coords(5,3)],'Color','black')
                        line([coords(2,1) coords(6,1)],[coords(2,2) coords(6,2)],[coords(2,3) coords(6,3)],'Color','black')
                        line([coords(3,1) coords(7,1)],[coords(3,2) coords(7,2)],[coords(3,3) coords(7,3)],'Color','black')
                        line([coords(4,1) coords(8,1)],[coords(4,2) coords(8,2)],[coords(4,3) coords(8,3)],'Color','black')


                    otherwise
                        error('not made yet')
                end
            end
        end
        %% Mesh Writer
        function MeshWrite(obj,gmsh_filename)
            %% MeshClass.MeshWrite
            %  Write external GMSH mesh files
            %
            % Syntax:
            %   obj = obj.MeshWrite( gmsh_filename );
            %
            % Inputs:
            %   gmsh_filename : filename
            %
            % Outputs:
            %
            % See also,
            gmsh = fopen ( gmsh_filename, 'wt' );

            if ( gmsh < 0 )
                fprintf ( 1, '\n' );
                fprintf ( 1, 'GMSH_MESH3D_WRITE - Error!\n' );
                fprintf ( 1, '  Could not open the output file.\n' );
                error ( 'GMSH_MESH3D_WRITE - Error!' );
            end
            %
            %  Write the data.
            %
            fprintf ( gmsh, '$MeshFormat\n' );
            fprintf ( gmsh, '2.2 0 8\n' );
            fprintf ( gmsh, '$EndMeshFormat\n' );

            fprintf ( gmsh, '$Nodes\n' );
            fprintf ( gmsh, '%d\n', size(obj.Coords,1) );
            for node = 1 : size(obj.Coords,1)
                fprintf ( gmsh, '%d %g %g %g\n', node,obj.Coords(node,1),obj.Coords(node,2),obj.Coords(node,3));
            end
            fprintf ( gmsh, '$EndNodes\n' );
            %
            %  These are the Gmsh codes for 4, 10 and 20 node tetrahedral elements.
            %
            %             if ( element_order == 4 )
            %                 element_type = 4;
            %             elseif ( element_order == 10 )
            %                 element_type = 11;
            %             elseif ( element_order == 20 )
            %                 element_type = 29;
            %             end
            element_type = 5;

            tag_num = 2;
            tag1 = 0;
            fprintf ( gmsh, '$Elements\n' );
            fprintf ( gmsh, '%d\n', size(obj.Connectivity,1) );
            for element = 1 : size(obj.Connectivity,1)
                fprintf ( gmsh, '%d  %d  %d  %d  %d', ...
                    element, element_type, tag_num, tag1, element );
                for vertex = 1 : 8
                    fprintf ( gmsh, '  %d', obj.Connectivity(element,vertex+1) );
                end
                fprintf ( gmsh, '\n' );
            end
            fprintf ( gmsh, '$EndElements\n' );

            fclose ( gmsh );
        end
    end
    methods(Static)
        function order = RotateHex(rot, sym)
            %% MeshClass.RotateHex
            %   Rotate the hexa to get some order in the template hexa ref
            %
            % Syntax:
            %   MeshClass.RotateHex (  )
            %
            % Inputs:
            %
            % Outputs:
            %
            % See also,

            order = 1:8;
            if rot(1)==-1
                order = order([3,2,6,7,0,1,5,4]+1);
            elseif rot(1)==1
                order = order([4,5,1,0,7,6,2,3]+1);
            end
            if rot(2)==-1
                order = order([4,0,3,7,5,1,2,6]+1);
            elseif rot(2)==1
                order = order([1,5,6,2,0,4,7,3]+1);
            end
            if rot(3)==-1
                order = order([1,2,3,0,5,6,7,4]+1);
            elseif rot(3)==1
                order = order([3,0,1,2,7,4,5,6]+1);
            end


            if(sym(1)==1)
                order = order([1,0,3,2,5,4,7,6]+1);
            end
            if(sym(2)==1)
                order = order([3,2,1,0,7,6,5,4]+1);
            end
            if(sym(3)==1)
                order = order([4,5,6,7,0,1,2,3]+1);
            end
        end
    end
end