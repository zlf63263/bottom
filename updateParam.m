% Update variable if designated as a swept parameter
function updatedValue = updateParam(origValue, parA, parB, parC)

% Keep original value if not a parameter
if( isnumeric(origValue) )
    updatedValue = origValue;
% Check if variable is parameter A
elseif( strcmp(origValue, 'param_A') )
    updatedValue = parA;
% Check if variable is parameter B
elseif( strcmp(origValue,'param_B') )
    updatedValue = parB;  
% Check if variable is parameter C
elseif( strcmp(origValue,'param_C') )
    updatedValue = parC;
end               

end

