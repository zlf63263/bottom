% Load imaginary component from library
function imaginary = imaginaryComponent(layer)

% Load library of materials
library = load('library-Diane.mat');

% Set layer to null if not present
if( strcmp(layer, 'N/A') )
    imaginary = double.empty(451,0);
% Set layer to corresponding material
else
    imaginary = library.(strcat(layer, '_k'));
end           

end
