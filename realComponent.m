% Load real component from library
function real = realComponent(layer)

% Load library of materials
library = load('library-Diane.mat');

% Set layer to null if not present
if( strcmp(layer, 'N/A') )
    real = double.empty(451,0);
% Set layer to corresponding material
else
    real = library.(strcat(layer, '_n'));
end           

end

