clear all; close all; clc;

% Initialize storage
nodeCoords = [];
elemList = [];
elemCounter = 1;

filesAndNodes = {
    'Beams/RedBeam.csv', 41, 1;
    'Beams/GreenBeams.csv', 6, 2;
    'Beams/BlueBeam1.csv', 8, 3;
    'Beams/BlueBeam2.csv', 6, 3;
    'Beams/BlueBeam3.csv', 8, 3;
};

% Tolerance for node uniqueness
tolerance = 1e-4;  % Adjust tolerance to avoid duplicate nodes
roundTol = 1e-4;   % Tolerance for deduplication

for f = 1:size(filesAndNodes, 1)
    fileNameWithPath = filesAndNodes{f, 1};
    nodesPerElement = filesAndNodes{f, 2};
    elemProperty = filesAndNodes{f, 3};
    
    warning('off', 'all');
    T = readtable(fileNameWithPath);
    warning('on', 'all');

    for i = 1:height(T)
        p1 = [T.StartX(i), T.StartY(i)];
        p2 = [T.EndX(i), T.EndY(i)];

        lineNodeIDs = zeros(nodesPerElement, 1);

        for j = 0:(nodesPerElement - 1)
            alpha = j / (nodesPerElement - 1);
            pt = (1 - alpha) * p1 + alpha * p2;

            % Check distances to all existing nodes using pdist2
            if isempty(nodeCoords)
                dists = [];
            else
                dists = pdist2(nodeCoords, pt);
            end

            matchIdx = find(dists < tolerance, 1);  % Find if there's any close enough match

            if ~isempty(matchIdx)  % If a match is found, use the existing node ID
                nodeID = matchIdx;
            else
                nodeCoords(end+1, :) = pt;  % Add new node if no match
                nodeID = size(nodeCoords, 1);
            end

            lineNodeIDs(j+1) = nodeID;
        end

        for j = 1:(nodesPerElement - 1)
            elemList(end+1,:) = [elemCounter, lineNodeIDs(j), lineNodeIDs(j+1), elemProperty];
            elemCounter = elemCounter + 1;
        end
    end
end

% --- Visualize the node positions ---
figure;
scatter(nodeCoords(:,1), nodeCoords(:,2), 'filled');
title('Node Positions');
xlabel('X Coordinate (m)');
ylabel('Y Coordinate (m)');
grid on;

% --- Deduplicate nodeCoords using a tolerance ---
roundedCoords = round(nodeCoords / roundTol) * roundTol;  % Round coordinates for precision

% Find unique coordinates and index mapping
[uniqueCoords, ~, idxMap] = unique(roundedCoords, 'rows', 'stable');

% Update element connectivity to use new unique node indices
elemList(:,2) = idxMap(elemList(:,2));
elemList(:,3) = idxMap(elemList(:,3));

% --- Write unique nodes to Excel ---
if isfile('inputs.xlsx')
    delete('inputs.xlsx');
end

nodeIDs = (1:size(uniqueCoords,1))';
constr_code = zeros(size(uniqueCoords,1), 1);
nodeTable = table(nodeIDs, constr_code(:), constr_code(:), constr_code(:),...
        uniqueCoords(:,1), uniqueCoords(:,2), ...
        'VariableNames', {'NodeID','X_constr','Y_constr','T_constr', 'X', 'Y'});
writetable(nodeTable, 'inputs.xlsx', 'Sheet', 'Nodes');

% --- Write updated elements to Excel ---
elementIDs = (1:size(elemList,1))';
elementTable = table(elementIDs, elemList(:,2), elemList(:,3), elemList(:,4), ...
    'VariableNames', {'ElementID', 'StartNode', 'EndNode', 'Property'});
writetable(elementTable, 'inputs.xlsx', 'Sheet', 'Elements');

% --- Optional: Sanity check for any residual duplicates ---
D = pdist2(uniqueCoords, uniqueCoords);
D(logical(eye(size(D)))) = NaN;  % Ignore self-distances
[minDist, idx] = min(D(:));

if minDist < roundTol
    [i, j] = ind2sub(size(D), idx);
    fprintf("❌ Warning: Nodes %d and %d are %.10f units apart (too close)\n", i, j, minDist);
else
    fprintf("✅ Export complete. %d unique nodes written. No close duplicates found.\n", size(uniqueCoords,1));
end
